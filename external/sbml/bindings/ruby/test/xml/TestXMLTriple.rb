#
# @file    TestXMLTriple.rb
# @brief   XMLTriple unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Michael Hucka <mhucka@caltech.edu> 
#
# $Id: TestXMLTriple.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/xml/TestXMLTriple.rb $
#
# This test file was converted from src/sbml/test/TestXMLTriple.c
# with the help of conversion sciprt (ctest_converter.pl).
#
#<!---------------------------------------------------------------------------
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright 2005-2009 California Institute of Technology.
# Copyright 2002-2005 California Institute of Technology and
#                     Japan Science and Technology Corporation.
# 
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
#--------------------------------------------------------------------------->*/
require 'test/unit'
require 'libSBML'

class TestXMLTriple < Test::Unit::TestCase

  def test_XMLTriple_create
    t = LibSBML::XMLTriple.new()
    assert( t != nil )
    assert( t.isEmpty() != false )
    t = nil
    t = LibSBML::XMLTriple.new("attr", "uri", "prefix")
    assert( (  "attr" != t.getName() ) == false )
    assert( (  "uri" != t.getURI() ) == false )
    assert( (  "prefix" != t.getPrefix() ) == false )
    assert( (  "prefix:attr" != t.getPrefixedName() ) == false )
    assert( t.isEmpty() == false )
    t = nil
    t = LibSBML::XMLTriple.new("attr", "uri", "")
    assert( (  "attr" != t.getName() ) == false )
    assert( (  "uri" != t.getURI() ) == false )
    assert( t.getPrefix() == "" )
    assert( (  "attr" != t.getPrefixedName() ) == false )
    assert( t.isEmpty() == false )
    t = nil
  end

end
