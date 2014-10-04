$:.unshift File.expand_path("../lib", __FILE__)
require 'rake'
require 'biotcm-cspn'

Gem::Specification.new do |s|
  s.name        = "biotcm-cspn"
  s.platform    = Gem::Platform::RUBY
  s.summary     = "A method to build pathway interaction network"
  s.description = "A method to build pathway interaction network."

  s.version     = BioTCM::Apps::CSPN::VERSION
  s.license     = 'MIT'

  s.authors     = ["Peng Zhang", "Aidi Stan"]
  s.email       = ["zp14@mails.tsinghua.edu.cn", "aidistan@live.cn"]
  s.homepage    = "http://biotcm.github.io/biotcm-cspn/"

  s.files         = FileList['lib/**/*', 'ref/**/*'].to_a
  s.require_paths = ['lib']

  s.required_ruby_version = '>= 2.0.0'
  s.add_runtime_dependency "biotcm", "~> 0.1.0"  
end
