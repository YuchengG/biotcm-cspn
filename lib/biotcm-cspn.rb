require 'biotcm'
require 'optparse'

# CSPN
class BioTCM::Apps::CSPN < BioTCM::Apps::App
  include BioTCM::Modules::WorkingDir

  # Version
  VERSION = "0.0.1"
  # Permutation times
  PERMUTATION = 1000
  # Exchange times
  EXCHANGE = 100000
  # Alpha level
  ALPHA = 0.05

  # Set wd
  def initialize
    self.wd = File.expand_path('../..', __FILE__)
  end
  # Run
  def run

    ###
    # Get options
    #

    options = {
      mode:'or',
      ppi:path_to('ref/ppi/combined.txt'),
      pathway_list:path_to('ref/pathway/default.txt'),
    }
    OptionParser.new do |opts|
      opts.banner = "Usage: biotcm cspn [OPTIONS]"

      opts.on("--gene-list [FILE]", String, "a list of genes") do |v|
        options[:gene_list] = v
      end
      opts.on("--ppi [FILE]", String, "a PPI file") do |v|
        options[:ppi] = v
      end
      opts.on("--pathway-list [FILE]", String, "a list of concerned KEGG pathways") do |v|
        options[:pathway_list] = v
      end
      opts.on("--mode [MODE]", String, "and, or") do |v|
        options[:mode] = 'and' if v =~ /and/i
        options[:mode] = 'or'  if v =~ /or/i
      end
    end.parse!
    raise 'Please give a list of genes at least' unless options[:gene_list]

    ###
    # Preparation
    #
    
    # Active PPI
    active_ppi = {}
    begin
      # Gene list
      gene = {}
      File.open(options[:gene_list]).each { |l| gene[l.chomp] = true }
      # Filter ppi
      File.open(options[:ppi]).each do |line|
        col = line.split("\t")
        if options[:mode] == 'and'
          next unless gene[col[0]] && gene[col[1]]
        else
          next unless gene[col[0]] || gene[col[1]]
        end
        active_ppi[col[0..1].join(' ')] = true
      end
    end

    # Genes of all pathways
    pathway = {}
    Dir[path_to "ref/pathway/gene_list/*.txt"].each do |filepath|
      /(?<id>hsa\d+)/ =~ filepath
      pathway[id] = File.open(filepath).collect { |l| l.chomp }
    end

    # Pathway list
    pathways = File.open(options[:pathway_list]).collect { |l| l.chomp }
    if (missing_pathways = pathways - pathway.keys).size > 0
      BioTCM.logger.warn('CSPN') do
        'Ignore unknown pathways: ' + missing_pathways.join(', ')
      end
      pathways -= missing_pathways
    end

    # Make result dir
    Dir.mkdir('cspn') unless Dir.exists?('cspn')
    Dir.mkdir('cspn/temp') unless Dir.exists?('cspn/temp')
    Dir.mkdir('cspn/result') unless Dir.exists?('cspn/result')

    ###
    # Significantlly interacted pathways
    #

    # Hash to store counts
    pathway_pair = {}
    pathways.product(pathways) do |p_p|
      next if p_p[0] == p_p[1]
      next if pathway_pair[p_p.reverse]
      pathway_pair[p_p] = []
    end

    # Real counts
    pathway_pair.each_key do |p_p|
      count = 0
      pathway[p_p[0]].product(pathway[p_p[1]]).each do |p_g|
        count += 1 if active_ppi[p_g.join(' ')] || active_ppi[p_g.reverse.join(' ')]
      end
      pathway_pair[p_p]<<count
    end
    File.open('cspn/temp/pathways_0.txt', 'w').puts pathways.collect { |w| [w, pathway[w]].join("\t") }
    File.open('cspn/temp/counts_0.txt', 'w').puts pathway_pair.to_a.collect { |a| a.join("\t") }

    # Permutation
    _pathway = pathway # Stash
    PERMUTATION.times do |permutation_time|
      BioTCM.logger.info('CSPN') { "Permutation #{permutation_time + 1}" }
      pathway = {}
      pathways.each { |path| pathway[path] = _pathway[path].clone }

      EXCHANGE.times do
        p1, p2 = pathways.shuffle[0..1]
        s1 = rand(pathway[p1].size)
        s2 = rand(pathway[p2].size) 
        while (pathway[p1][s1] == pathway[p2][s2])
          s1 = rand(pathway[p1].size)
          s2 = rand(pathway[p2].size)
        end
        pathway[p1][s1], pathway[p2][s2] = pathway[p2][s2], pathway[p1][s1]
      end

      pathway_pair.each_key do |p_p|
        count = 0
        pathway[p_p[0]].product(pathway[p_p[1]]).each do |p_g|
          count += 1 if active_ppi[p_g.join(' ')] || active_ppi[p_g.reverse.join(' ')]
        end
        pathway_pair[p_p]<<count
      end
      File.open("cspn/temp/pathways_#{permutation_time + 1}.txt", 'w').puts pathway.to_a.collect { |a| a.join("\t") }
      File.open("cspn/temp/counts_#{permutation_time + 1}.txt", 'w') do |fout|
        pathway_pair.each { |k, v| fout.puts [k, v.last].join("\t") }
      end
    end
    pathway = _pathway # Pop
    
    # Calculate empirical p-values
    BioTCM.logger.info('CSPN') { 'Calculating empirical p-values' }
    pathway_pair.each do |k, v|
      pathway_pair[k] = 1 - (v.sort.index(v[0]) + (v.count(v[0])-1.0) /2) / (PERMUTATION + 1)
    end
    File.open('cspn/temp/p-values.txt', 'w').puts pathway_pair.to_a.collect { |a| a.join("\t") }

    # Adjust p-values
    `Rscript #{path_to 'lib/biotcm-cspn/adjust_p_value.r'}`
    sleep(3) # Make sure the file has already be flushed

    ###
    # Output results
    #

    # Pathway-pathway interaction
    ppis_pathway = File.open('cspn/temp/adj.p-values.txt').collect do |line|
      col = line.chomp.split("\t")
    end.select { |a| a[2].to_f < ALPHA }

    File.open('cspn/result/pathway_pathway_interactions.txt', 'w') do |fout|
      fout.puts ppis_pathway.collect { |a| a.join("\t") }
    end

    # Active PPIs
    ppis_gene = {}
    ppis_pathway.each do |p_p|
      str_p_p = p_p[0..1].join('-')
      pathway[p_p[0]].product(pathway[p_p[1]]).each do |p_g|
        if active_ppi[p_g.join(' ')]
          ppis_gene[p_g.join(' ')] ? ppis_gene[p_g.join(' ')]<<str_p_p : ppis_gene[p_g.join(' ')] = [str_p_p]
        elsif active_ppi[p_g.reverse.join(' ')]
          ppis_gene[p_g.reverse.join(' ')] ? ppis_gene[p_g.reverse.join(' ')]<<str_p_p : ppis_gene[p_g.reverse.join(' ')] = [str_p_p]
        end
      end
    end
    File.open('cspn/result/active_ppi_involved.txt', 'w') do |fout|
      ppis_gene.keys.sort_by { |k| -ppis_gene[k].size }.each do |p_g|
        fout.puts p_g.split(' ').push(ppis_gene[p_g].join(', ')).join("\t")
      end
    end
  end
end
