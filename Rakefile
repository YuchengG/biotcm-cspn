$:.push File.expand_path("../lib", __FILE__)
require 'biotcm-cspn'

# clean
desc "Clean the directory"
task :clean do
  (FileList[".yardoc", "doc", "*.gem", "test/cspn"]).each do |d|
    FileUtils.rm_r(d) rescue nil
  end
end

# gem
desc "Build the gem"
task :gem do
  system("gem build #{File.dirname(__FILE__)}" + "/biotcm-cspn.gemspec")
end

# install
desc "Install the gem"
task :install => :gem do
  system("gem install biotcm-cspn-#{BioTCM::Apps::CSPN::VERSION}.gem --local --no-document")
end

# uninstall
desc "Uninstall the gem"
task :uninstall do
  system("gem uninstall biotcm-cspn")
end
