#
# @file    TestSBMLDocument.rb
# @brief   SBMLDocument unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Ben Bornstein 
#
# $Id: TestSBMLDocument.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestSBMLDocument.rb $
#
# This test file was converted from src/sbml/test/TestSBMLDocument.c
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

class TestSBMLDocument < Test::Unit::TestCase

  def test_SBMLDocument_create
    d = LibSBML::SBMLDocument.new()
    assert( d.getTypeCode() == LibSBML::SBML_DOCUMENT )
    assert( d.getNotes() == nil )
    assert( d.getAnnotation() == nil )
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 4 )
    assert( d.getNumErrors() == 0 )
    d = nil
  end

  def test_SBMLDocument_createWith
    d = LibSBML::SBMLDocument.new(1,2)
    assert( d.getTypeCode() == LibSBML::SBML_DOCUMENT )
    assert( d.getNotes() == nil )
    assert( d.getAnnotation() == nil )
    assert( d.getLevel() == 1 )
    assert( d.getVersion() == 2 )
    assert( d.getNumErrors() == 0 )
    d = nil
  end

  def test_SBMLDocument_free_NULL
  end

  def test_SBMLDocument_setLevelAndVersion
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,2,false)
    m1 = LibSBML::Model.new(2,2)
    d.setModel(m1)
    assert( d.setLevelAndVersion(2,3,false) == true )
    assert( d.setLevelAndVersion(2,1,false) == true )
    assert( d.setLevelAndVersion(1,2,false) == true )
    assert( d.setLevelAndVersion(1,1,false) == false )
    d = nil
  end

  def test_SBMLDocument_setLevelAndVersion_Error
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,1,false)
    m1 = LibSBML::Model.new(2,1)
    u = LibSBML::Unit.new(2,1)
    u.setKind(LibSBML::UnitKind_forName("mole"))
    u.setOffset(3.2)
    ud = LibSBML::UnitDefinition.new(2,1)
    ud.setId( "ud")
    ud.addUnit(u)
    m1.addUnitDefinition(ud)
    d.setModel(m1)
    assert( d.setLevelAndVersion(2,2,false) == false )
    assert( d.setLevelAndVersion(2,3,false) == false )
    assert( d.setLevelAndVersion(1,2,false) == false )
    assert( d.setLevelAndVersion(1,1,false) == false )
    d = nil
  end

  def test_SBMLDocument_setLevelAndVersion_UnitsError
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,4,false)
    m1 = d.createModel()
    c = m1.createCompartment()
    c.setId( "c")
    p = m1.createParameter()
    p.setId( "p")
    p.setUnits( "mole")
    r = m1.createAssignmentRule()
    r.setVariable( "c")
    r.setFormula( "p*p")
    assert( d.setLevelAndVersion(2,2,false) == true )
    assert( d.setLevelAndVersion(2,3,false) == true )
    assert( d.setLevelAndVersion(1,2,false) == true )
    assert( d.setLevelAndVersion(1,1,false) == false )
    d = nil
  end

  def test_SBMLDocument_setLevelAndVersion_Warning
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,2,false)
    m1 = LibSBML::Model.new(2,2)
    (m1).setSBOTerm(2)
    d.setModel(m1)
    assert( d.setLevelAndVersion(2,3,false) == true )
    assert( d.setLevelAndVersion(2,1,false) == true )
    assert( d.setLevelAndVersion(1,2,false) == true )
    assert( d.setLevelAndVersion(1,1,false) == false )
    d = nil
  end

  def test_SBMLDocument_setModel
    d = LibSBML::SBMLDocument.new()
    m1 = LibSBML::Model.new(2,4)
    m2 = LibSBML::Model.new(2,4)
    assert( d.getModel() == nil )
    d.setModel(m1)
    assert( d.getModel() != m1 )
    d.setModel(d.getModel())
    assert( d.getModel() != m1 )
    d.setModel(m2)
    assert( d.getModel() != m2 )
    d = nil
  end

  def test_SBMLDocument_setModel1
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,2,false)
    m1 = LibSBML::Model.new(2,1)
    i = d.setModel(m1)
    assert( i == LibSBML::LIBSBML_VERSION_MISMATCH )
    assert( d.getModel() == nil )
    d = nil
  end

  def test_SBMLDocument_setModel2
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,2,false)
    m1 = LibSBML::Model.new(1,2)
    i = d.setModel(m1)
    assert( i == LibSBML::LIBSBML_LEVEL_MISMATCH )
    assert( d.getModel() == nil )
    d = nil
  end

  def test_SBMLDocument_setModel3
    d = LibSBML::SBMLDocument.new()
    d.setLevelAndVersion(2,2,false)
    m1 = LibSBML::Model.new(2,2)
    i = d.setModel(m1)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( d.getModel() != nil )
    d = nil
  end

end
