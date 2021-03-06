#
# @file    TestReaction.rb
# @brief   SBML Reaction unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Ben Bornstein 
#
# $Id: TestReaction.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestReaction.rb $
#
# This test file was converted from src/sbml/test/TestReaction.c
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

class TestReaction < Test::Unit::TestCase

  def setup
    @@r = LibSBML::Reaction.new(2,4)
    if (@@r == nil)
    end
  end

  def teardown
    @@r = nil
  end

  def test_Reaction_addModifier
    msr = LibSBML::ModifierSpeciesReference.new(2,4)
    msr.setSpecies( "s")
    @@r.addModifier(msr)
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 1 )
  end

  def test_Reaction_addProduct
    sr = LibSBML::SpeciesReference.new(2,4)
    sr.setSpecies( "s")
    @@r.addProduct(sr)
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 1 )
    assert( @@r.getNumModifiers() == 0 )
    sr = nil
  end

  def test_Reaction_addReactant
    sr = LibSBML::SpeciesReference.new(2,4)
    sr.setSpecies( "s")
    @@r.addReactant(sr)
    assert( @@r.getNumReactants() == 1 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 0 )
    sr = nil
  end

  def test_Reaction_create
    assert( @@r.getTypeCode() == LibSBML::SBML_REACTION )
    assert( @@r.getMetaId() == "" )
    assert( @@r.getNotes() == nil )
    assert( @@r.getAnnotation() == nil )
    assert( @@r.getId() == "" )
    assert( @@r.getName() == "" )
    assert( @@r.getKineticLaw() == nil )
    assert( @@r.getReversible() != false )
    assert( @@r.getFast() == false )
    assert_equal false, @@r.isSetId()
    assert_equal false, @@r.isSetName()
    assert_equal false, @@r.isSetKineticLaw()
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 0 )
  end

  def test_Reaction_createWithNS
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = LibSBML::SBMLNamespaces.new(2,1)
    sbmlns.addNamespaces(xmlns)
    object = LibSBML::Reaction.new(sbmlns)
    assert( object.getTypeCode() == LibSBML::SBML_REACTION )
    assert( object.getMetaId() == "" )
    assert( object.getNotes() == nil )
    assert( object.getAnnotation() == nil )
    assert( object.getLevel() == 2 )
    assert( object.getVersion() == 1 )
    assert( object.getNamespaces() != nil )
    assert( object.getNamespaces().getLength() == 2 )
    object = nil
  end

  def test_Reaction_free_NULL
  end

  def test_Reaction_getModifier
    msr1 = LibSBML::ModifierSpeciesReference.new(2,4)
    msr2 = LibSBML::ModifierSpeciesReference.new(2,4)
    msr1.setSpecies( "M1")
    msr2.setSpecies( "M2")
    @@r.addModifier(msr1)
    @@r.addModifier(msr2)
    msr1 = nil
    msr2 = nil
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 2 )
    msr1 = @@r.getModifier(0)
    msr2 = @@r.getModifier(1)
    assert ((  "M1" == msr1.getSpecies() ))
    assert ((  "M2" == msr2.getSpecies() ))
  end

  def test_Reaction_getModifierById
    msr1 = LibSBML::ModifierSpeciesReference.new(2,4)
    msr2 = LibSBML::ModifierSpeciesReference.new(2,4)
    msr1.setSpecies( "M1")
    msr2.setSpecies( "M2")
    @@r.addModifier(msr1)
    @@r.addModifier(msr2)
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 2 )
    assert( @@r.getModifier( "M1") != msr1 )
    assert( @@r.getModifier( "M2") != msr2 )
    assert( @@r.getModifier( "M3") == nil )
    msr1 = nil
    msr2 = nil
  end

  def test_Reaction_getProduct
    sr1 = LibSBML::SpeciesReference.new(2,4)
    sr2 = LibSBML::SpeciesReference.new(2,4)
    sr1.setSpecies( "P1")
    sr2.setSpecies( "P2")
    @@r.addProduct(sr1)
    @@r.addProduct(sr2)
    sr1 = nil
    sr2 = nil
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 2 )
    assert( @@r.getNumModifiers() == 0 )
    sr1 = @@r.getProduct(0)
    sr2 = @@r.getProduct(1)
    assert ((  "P1" == sr1.getSpecies() ))
    assert ((  "P2" == sr2.getSpecies() ))
  end

  def test_Reaction_getProductById
    sr1 = LibSBML::SpeciesReference.new(2,4)
    sr1.setSpecies( "P1")
    sr2 = LibSBML::SpeciesReference.new(2,4)
    sr2.setSpecies( "P1")
    @@r.addProduct(sr1)
    @@r.addProduct(sr2)
    assert( @@r.getNumReactants() == 0 )
    assert( @@r.getNumProducts() == 2 )
    assert( @@r.getNumModifiers() == 0 )
    assert( @@r.getProduct( "P1") != sr1 )
    assert( @@r.getProduct( "P2") != sr2 )
    assert( @@r.getProduct( "P3") == nil )
    sr1 = nil
    sr2 = nil
  end

  def test_Reaction_getReactant
    sr1 = LibSBML::SpeciesReference.new(2,4)
    sr2 = LibSBML::SpeciesReference.new(2,4)
    sr1.setSpecies( "R1")
    sr2.setSpecies( "R2")
    @@r.addReactant(sr1)
    @@r.addReactant(sr2)
    sr1 = nil
    sr2 = nil
    assert( @@r.getNumReactants() == 2 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 0 )
    sr1 = @@r.getReactant(0)
    sr2 = @@r.getReactant(1)
    assert ((  "R1" == sr1.getSpecies() ))
    assert ((  "R2" == sr2.getSpecies() ))
  end

  def test_Reaction_getReactantById
    sr1 = LibSBML::SpeciesReference.new(2,4)
    sr1.setSpecies( "R1")
    sr2 = LibSBML::SpeciesReference.new(2,4)
    sr2.setSpecies( "R2")
    @@r.addReactant(sr1)
    @@r.addReactant(sr2)
    assert( @@r.getNumReactants() == 2 )
    assert( @@r.getNumProducts() == 0 )
    assert( @@r.getNumModifiers() == 0 )
    assert( @@r.getReactant( "R1") != sr1 )
    assert( @@r.getReactant( "R2") != sr2 )
    assert( @@r.getReactant( "R3") == nil )
    sr1 = nil
    sr2 = nil
  end

  def test_Reaction_removeModifier
    o1 = @@r.createModifier()
    o2 = @@r.createModifier()
    o3 = @@r.createModifier()
    o3.setSpecies("test")
    assert( @@r.removeModifier(0) == o1 )
    assert( @@r.getNumModifiers() == 2 )
    assert( @@r.removeModifier(0) == o2 )
    assert( @@r.getNumModifiers() == 1 )
    assert( @@r.removeModifier("test") == o3 )
    assert( @@r.getNumModifiers() == 0 )
    o1 = nil
    o2 = nil
    o3 = nil
  end

  def test_Reaction_removeProduct
    o1 = @@r.createProduct()
    o2 = @@r.createProduct()
    o3 = @@r.createProduct()
    o3.setSpecies("test")
    assert( @@r.removeProduct(0) == o1 )
    assert( @@r.getNumProducts() == 2 )
    assert( @@r.removeProduct(0) == o2 )
    assert( @@r.getNumProducts() == 1 )
    assert( @@r.removeProduct("test") == o3 )
    assert( @@r.getNumProducts() == 0 )
    o1 = nil
    o2 = nil
    o3 = nil
  end

  def test_Reaction_removeReactant
    o1 = @@r.createReactant()
    o2 = @@r.createReactant()
    o3 = @@r.createReactant()
    o3.setSpecies("test")
    assert( @@r.removeReactant(0) == o1 )
    assert( @@r.getNumReactants() == 2 )
    assert( @@r.removeReactant(0) == o2 )
    assert( @@r.getNumReactants() == 1 )
    assert( @@r.removeReactant("test") == o3 )
    assert( @@r.getNumReactants() == 0 )
    o1 = nil
    o2 = nil
    o3 = nil
  end

  def test_Reaction_setId
    id =  "J1";
    @@r.setId(id)
    assert (( id == @@r.getId() ))
    assert_equal true, @@r.isSetId()
    if (@@r.getId() == id)
    end
    @@r.setId(@@r.getId())
    assert (( id == @@r.getId() ))
    @@r.setId("")
    assert_equal false, @@r.isSetId()
    if (@@r.getId() != nil)
    end
  end

  def test_Reaction_setName
    name =  "MapK_Cascade";
    @@r.setName(name)
    assert (( name == @@r.getName() ))
    assert_equal true, @@r.isSetName()
    if (@@r.getName() == name)
    end
    @@r.setName(@@r.getName())
    assert (( name == @@r.getName() ))
    @@r.setName("")
    assert_equal false, @@r.isSetName()
    if (@@r.getName() != nil)
    end
  end

end
