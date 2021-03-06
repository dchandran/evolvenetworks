#
# @file    TestInitialAssignment.rb
# @brief   SBML InitialAssignment unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Sarah Keating 
#
# $Id: TestInitialAssignment.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestInitialAssignment.rb $
#
# This test file was converted from src/sbml/test/TestInitialAssignment.c
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

class TestInitialAssignment < Test::Unit::TestCase

  def setup
    @@ia = LibSBML::InitialAssignment.new(2,4)
    if (@@ia == nil)
    end
  end

  def teardown
    @@ia = nil
  end

  def test_InitialAssignment_create
    assert( @@ia.getTypeCode() == LibSBML::SBML_INITIAL_ASSIGNMENT )
    assert( @@ia.getMetaId() == "" )
    assert( @@ia.getNotes() == nil )
    assert( @@ia.getAnnotation() == nil )
    assert( @@ia.getSymbol() == "" )
    assert( @@ia.getMath() == nil )
  end

  def test_InitialAssignment_createWithNS
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = LibSBML::SBMLNamespaces.new(2,3)
    sbmlns.addNamespaces(xmlns)
    object = LibSBML::InitialAssignment.new(sbmlns)
    assert( object.getTypeCode() == LibSBML::SBML_INITIAL_ASSIGNMENT )
    assert( object.getMetaId() == "" )
    assert( object.getNotes() == nil )
    assert( object.getAnnotation() == nil )
    assert( object.getLevel() == 2 )
    assert( object.getVersion() == 3 )
    assert( object.getNamespaces() != nil )
    assert( object.getNamespaces().getLength() == 2 )
    object = nil
  end

  def test_InitialAssignment_free_NULL
  end

  def test_InitialAssignment_setMath
    math = LibSBML::parseFormula("2 * k")
    @@ia.setMath(math)
    math1 = @@ia.getMath()
    assert( math1 != nil )
    formula = LibSBML::formulaToString(math1)
    assert( formula != nil )
    assert ((  "2 * k" == formula ))
    assert( @@ia.getMath() != math )
    assert_equal true, @@ia.isSetMath()
    @@ia.setMath(@@ia.getMath())
    math1 = @@ia.getMath()
    assert( math1 != nil )
    formula = LibSBML::formulaToString(math1)
    assert( formula != nil )
    assert ((  "2 * k" == formula ))
    assert( @@ia.getMath() != math )
    @@ia.setMath(nil)
    assert_equal false, @@ia.isSetMath()
    if (@@ia.getMath() != nil)
    end
    math = nil
  end

  def test_InitialAssignment_setSymbol
    symbol =  "k2";
    @@ia.setSymbol(symbol)
    assert (( symbol == @@ia.getSymbol() ))
    assert_equal true, @@ia.isSetSymbol()
    if (@@ia.getSymbol() == symbol)
    end
    @@ia.setSymbol(@@ia.getSymbol())
    assert (( symbol == @@ia.getSymbol() ))
    @@ia.setSymbol("")
    assert_equal false, @@ia.isSetSymbol()
    if (@@ia.getSymbol() != nil)
    end
  end

end
