require 'biotcm'
require 'optparse'

# CSPN: A method to build pathway interaction network.
#
# == Reference
# Yezhou Huang and Shao Li, “Detection of characteristic sub pathway 
# network for angio-genesis based on the comprehensive pathway network”, 
# BMC Bioinformatics, Volume 11 Supplement 1, 2010: Selected articles 
# from the Eighth Asia-Pacific Bioinformatics Conference, Bangalore, 
# India, January 2010 (APBC 2010).
#
class BioTCM::Apps::CSPN < BioTCM::Apps::App
  include BioTCM::Modules::WorkingDir

  # Version
  VERSION = "0.0.1"
  # Permutation times
  PERMUTATION = 1
  # Exchange times
  EXCHANGE = 100000

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
    pathways.product(pathways) do |pair|
      next if pair[0] == pair[1]
      next if pathway_pair[pair.reverse]
      pathway_pair[pair] = []
    end

    # Real counts
    pathway_pair.each_key do |pair|
      count = 0
      pathway[pair[0]].product(pathway[pair[1]]).each do |pair|
        count += 1 if active_ppi[pair.join(' ')] || active_ppi[pair.reverse.join(' ')]
      end
      pathway_pair[pair]<<count
    end
    File.open('cspn/temp/pathways_0.txt', 'w').puts pathways.collect { |w| [w, pathway[w]].join("\t") }
    File.open('cspn/temp/counts_0.txt', 'w').puts pathway_pair.to_a.collect { |a| a.join("\t") }

    # Permutation
    _pathway = pathway
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

      pathway_pair.each_key do |pair|
        count = 0
        pathway[pair[0]].product(pathway[pair[1]]).each do |pair|
          count += 1 if active_ppi[pair.join(' ')] || active_ppi[pair.reverse.join(' ')]
        end
        pathway_pair[pair]<<count
      end
      File.open("cspn/temp/pathways_#{permutation_time + 1}.txt", 'w').puts pathway.to_a.collect { |a| a.join("\t") }
      File.open("cspn/temp/counts_#{permutation_time + 1}.txt", 'w') do |fout|
        pathway_pair.each { |k, v| fout.puts [k, v.last].join("\t") }
      end
    end
    
    # Calculate empirical pvalues
    BioTCM.logger.info('CSPN') { 'Calculating empirical p-values' }
    pathway_pair.each do |k, v|
      pathway_pair[k] = 1 - (v.sort.index(v[0]) + (v.count(v[0])-1.0) /2) / (PERMUTATION + 1)
    end
    File.open('cspn/result/p-values.txt', 'w').puts pathway_pair.to_a.collect { |a| a.join("\t") }

    `Rscript #{path_to 'lib/adjust_p_value.r'}`    

    ###
    # Output results
    #

  end
end
