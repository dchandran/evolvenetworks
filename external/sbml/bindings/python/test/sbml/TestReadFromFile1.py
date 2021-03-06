#
# @file    TestReadFromFile1.py
# @brief   Reads tests/l1v1-branch.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestReadFromFile1.py 8704 2009-01-04 02:26:05Z mhucka $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestReadFromFile1.py $
#
# This test file was converted from src/sbml/test/TestReadFromFile1.c
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
import sys
import unittest
import libsbml

class TestReadFromFile1(unittest.TestCase):


  def test_read_l1v1_branch(self):
    filename = "../../sbml/test/test-data/l1v1-branch.xml"
    d = libsbml.readSBML(filename)
    if (d == None):
      pass    
    self.assert_( d.getLevel() == 1 )
    self.assert_( d.getVersion() == 1 )
    m = d.getModel()
    self.assert_((  "Branch" == m.getName() ))
    self.assert_( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    self.assert_((  "compartmentOne" == c.getName() ))
    self.assert_( c.getVolume() == 1 )
    ud = c.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 1 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_LITRE )
    self.assert_( m.getNumSpecies() == 4 )
    s = m.getSpecies(0)
    self.assert_((  "S1"              == s.getName() ))
    self.assert_((  "compartmentOne"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == False )
    ud = s.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 2 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_MOLE )
    self.assert_( ud.getUnit(0).getExponent() == 1 )
    self.assert_( ud.getUnit(1).getKind() == libsbml.UNIT_KIND_LITRE )
    self.assert_( ud.getUnit(1).getExponent() == -1 )
    s = m.getSpecies(1)
    self.assert_((  "X0"              == s.getName() ))
    self.assert_((  "compartmentOne"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == True )
    s = m.getSpecies(2)
    self.assert_((  "X1"              == s.getName() ))
    self.assert_((  "compartmentOne"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == True )
    s = m.getSpecies(3)
    self.assert_((  "X2"              == s.getName() ))
    self.assert_((  "compartmentOne"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == True )
    self.assert_( m.getNumReactions() == 3 )
    r = m.getReaction(0)
    self.assert_((  "reaction_1" == r.getName() ))
    self.assert_( r.getReversible() == False )
    self.assert_( r.getFast() == False )
    ud = r.getKineticLaw().getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 2 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_MOLE )
    self.assert_( ud.getUnit(0).getExponent() == 1 )
    self.assert_( ud.getUnit(1).getKind() == libsbml.UNIT_KIND_LITRE )
    self.assert_( ud.getUnit(1).getExponent() == -1 )
    self.assert_( r.getKineticLaw().containsUndeclaredUnits() == True )
    r = m.getReaction(1)
    self.assert_((  "reaction_2" == r.getName() ))
    self.assert_( r.getReversible() == False )
    self.assert_( r.getFast() == False )
    r = m.getReaction(2)
    self.assert_((  "reaction_3" == r.getName() ))
    self.assert_( r.getReversible() == False )
    self.assert_( r.getFast() == False )
    r = m.getReaction(0)
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "X0" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "S1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    kl = r.getKineticLaw()
    self.assert_((  "k1 * X0" == kl.getFormula() ))
    self.assert_( kl.getNumParameters() == 1 )
    p = kl.getParameter(0)
    self.assert_((  "k1" == p.getName() ))
    self.assert_( p.getValue() == 0 )
    r = m.getReaction(1)
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "S1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "X1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    kl = r.getKineticLaw()
    self.assert_((  "k2 * S1" == kl.getFormula() ))
    self.assert_( kl.getNumParameters() == 1 )
    p = kl.getParameter(0)
    self.assert_((  "k2" == p.getName() ))
    self.assert_( p.getValue() == 0 )
    r = m.getReaction(2)
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "S1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "X2" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    kl = r.getKineticLaw()
    self.assert_((  "k3 * S1" == kl.getFormula() ))
    self.assert_( kl.getNumParameters() == 1 )
    p = kl.getParameter(0)
    self.assert_((  "k3" == p.getName() ))
    self.assert_( p.getValue() == 0 )
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestReadFromFile1))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
