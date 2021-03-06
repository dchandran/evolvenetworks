#
# @file    TestSBMLConvert.rb
# @brief   SBMLConvert unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Ben Bornstein 
#
# $Id: TestSBMLConvert.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestSBMLConvert.rb $
#
# This test file was converted from src/sbml/test/TestSBMLConvert.c
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

class TestSBMLConvert < Test::Unit::TestCase

  def test_SBMLConvert_addModifiersToReaction
    d = LibSBML::SBMLDocument.new(1,2)
    m = d.createModel()
    r = m.createReaction()
    kl = r.createKineticLaw()
    kl.setFormula( "k1*S1*S2*S3*S4*S5")
    s1 = m.createSpecies()
    s1.setId( "S1" )
    s2 = m.createSpecies()
    s2.setId( "S2")
    s3 = m.createSpecies()
    s3.setId( "S3")
    s4 = m.createSpecies()
    s4.setId( "S4")
    s5 = m.createSpecies()
    s5.setId( "S5")
    sr1 = r.createReactant()
    sr2 = r.createReactant()
    sr3 = r.createProduct()
    sr1.setSpecies( "S1")
    sr2.setSpecies( "S2")
    sr3.setSpecies( "S5")
    assert( r.getNumModifiers() == 0 )
    d.setLevelAndVersion(2,1,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 1 )
    assert( m.getReaction(0).getNumModifiers() == 2 )
    ssr1 = m.getReaction(0).getModifier(0)
    ssr2 = m.getReaction(0).getModifier(1)
    assert ((  "S3" == ssr1.getSpecies() ))
    assert ((  "S4" == ssr2.getSpecies() ))
    d = nil
  end

  def test_SBMLConvert_convertToL1_SBMLDocument
    d = LibSBML::SBMLDocument.new(2,1)
    d.setLevelAndVersion(1,2,false)
    assert( d.getLevel() == 1 )
    assert( d.getVersion() == 2 )
    d = nil
  end

  def test_SBMLConvert_convertToL1_Species_Amount
    d = LibSBML::SBMLDocument.new(2,1)
    m = d.createModel()
    sid =  "C";
    c = LibSBML::Compartment.new(2,4)
    s = LibSBML::Species.new(2,4)
    c.setId(sid)
    m.addCompartment(c)
    s.setCompartment(sid)
    s.setInitialAmount(2.34)
    m.addSpecies(s)
    d.setLevelAndVersion(1,2,false)
    assert( s.getInitialAmount() == 2.34 )
    d = nil
  end

  def test_SBMLConvert_convertToL1_Species_Concentration
    d = LibSBML::SBMLDocument.new(2,1)
    m = d.createModel()
    sid =  "C";
    c = LibSBML::Compartment.new(2,1)
    s = LibSBML::Species.new(2,1)
    c.setId(sid)
    c.setSize(1.2)
    m.addCompartment(c)
    s.setId( "s"  )
    s.setCompartment(sid)
    s.setInitialConcentration(2.34)
    m.addSpecies(s)
    d.setLevelAndVersion(1,2,false)
    s1 = m.getSpecies(0)
    assert( s1 != nil )
    assert ((  "C" == s1.getCompartment() ))
    assert( m.getCompartment( "C").getSize() == 1.2 )
    assert( s1.getInitialConcentration() == 2.34 )
    assert( s1.isSetInitialConcentration() == true )
    d = nil
  end

  def test_SBMLConvert_convertToL2_SBMLDocument
    d = LibSBML::SBMLDocument.new(1,2)
    d.setLevelAndVersion(2,1,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 1 )
    d.setLevelAndVersion(2,2,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 2 )
    d.setLevelAndVersion(2,3,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 3 )
    d = nil
  end

  def test_SBMLConvert_convertToL2v4_DuplicateAnnotations_doc
    d = LibSBML::SBMLDocument.new(2,1)
    d.createModel()
    annotation =  "<rdf/>\n<rdf/>";
    i = (d).setAnnotation(annotation)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 1 )
    assert( (d).getAnnotation().getNumChildren() == 2 )
    d.setLevelAndVersion(2,4,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 4 )
    assert( (d).getAnnotation().getNumChildren() == 1 )
    d = nil
  end

  def test_SBMLConvert_convertToL2v4_DuplicateAnnotations_model
    d = LibSBML::SBMLDocument.new(2,1)
    m = d.createModel()
    annotation =  "<rdf/>\n<rdf/>";
    i = (m).setAnnotation(annotation)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 1 )
    assert( (m).getAnnotation().getNumChildren() == 2 )
    d.setLevelAndVersion(2,4,false)
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 4 )
    m = d.getModel()
    assert( (m).getAnnotation().getNumChildren() == 1 )
    d = nil
  end

end
