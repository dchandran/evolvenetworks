#
# @file    TestStoichiometryMath.rb
# @brief   SBML StoichiometryMath unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Sarah Keating 
#
# $Id: TestStoichiometryMath.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestStoichiometryMath.rb $
#
# This test file was converted from src/sbml/test/TestStoichiometryMath.c
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

class TestStoichiometryMath < Test::Unit::TestCase

  def setup
    @@d = LibSBML::StoichiometryMath.new(2,4)
    if (@@d == nil)
    end
  end

  def teardown
    @@d = nil
  end

  def test_StoichiometryMath_create
    assert( @@d.getTypeCode() == LibSBML::SBML_STOICHIOMETRY_MATH )
    assert( @@d.getMetaId() == "" )
    assert( @@d.getNotes() == nil )
    assert( @@d.getAnnotation() == nil )
    assert( @@d.getMath() == nil )
  end

  def test_StoichiometryMath_createWithNS
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = LibSBML::SBMLNamespaces.new(2,1)
    sbmlns.addNamespaces(xmlns)
    object = LibSBML::StoichiometryMath.new(sbmlns)
    assert( object.getTypeCode() == LibSBML::SBML_STOICHIOMETRY_MATH )
    assert( object.getMetaId() == "" )
    assert( object.getNotes() == nil )
    assert( object.getAnnotation() == nil )
    assert( object.getLevel() == 2 )
    assert( object.getVersion() == 1 )
    assert( object.getNamespaces() != nil )
    assert( object.getNamespaces().getLength() == 2 )
    object = nil
  end

  def test_StoichiometryMath_free_NULL
  end

  def test_StoichiometryMath_setMath
    math = LibSBML::parseFormula("lambda(x, x^3)")
    @@d.setMath(math)
    math1 = @@d.getMath()
    assert( math1 != nil )
    formula = LibSBML::formulaToString(math1)
    assert( formula != nil )
    assert ((  "lambda(x, x^3)" == formula ))
    assert( @@d.getMath() != math )
    assert_equal true, @@d.isSetMath()
    @@d.setMath(@@d.getMath())
    math1 = @@d.getMath()
    assert( math1 != nil )
    formula = LibSBML::formulaToString(math1)
    assert( formula != nil )
    assert ((  "lambda(x, x^3)" == formula ))
    @@d.setMath(nil)
    assert_equal false, @@d.isSetMath()
    if (@@d.getMath() != nil)
    end
  end

  def test_StoichiometryMath_setMath1
    math = LibSBML::parseFormula("2 * k")
    i = @@d.setMath(math)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@d.getMath() != math )
    assert_equal true, @@d.isSetMath()
    i = @@d.setMath(nil)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@d.getMath() == nil )
    assert_equal false, @@d.isSetMath()
    math = nil
  end

  def test_StoichiometryMath_setMath2
    math = LibSBML::ASTNode.new(LibSBML::AST_TIMES)
    i = @@d.setMath(math)
    assert( i == LibSBML::LIBSBML_INVALID_OBJECT )
    assert_equal false, @@d.isSetMath()
    math = nil
  end

end
