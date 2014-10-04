#!/usr/bin/env ruby
$:.unshift File.expand_path("../../lib", __FILE__)
require 'biotcm-cspn'

ARGV<<'--gene-list'<<'list_short.txt'
BioTCM::Apps::CSPN.new.run
