#
# @file    TestReadFromFile3.py
# @brief   Reads tests/l1v1-rules.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestReadFromFile3.py 8704 2009-01-04 02:26:05Z mhucka $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestReadFromFile3.py $
#
# This test file was converted from src/sbml/test/TestReadFromFile3.c
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

class TestReadFromFile3(unittest.TestCase):


  def test_read_l1v1_rules(self):
    filename = "../../sbml/test/test-data/l1v1-rules.xml"
    d = libsbml.readSBML(filename)
    if (d == None):
      pass    
    self.assert_( d.getLevel() == 1 )
    self.assert_( d.getVersion() == 1 )
    m = d.getModel()
    self.assert_( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    self.assert_((  "cell" == c.getName() ))
    self.assert_( c.getVolume() == 1 )
    self.assert_( m.getNumSpecies() == 6 )
    s = m.getSpecies(0)
    self.assert_((  "s1"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 4 )
    self.assert_( s.getBoundaryCondition() == False )
    s = m.getSpecies(1)
    self.assert_((  "s2"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 2 )
    self.assert_( s.getBoundaryCondition() == False )
    s = m.getSpecies(2)
    self.assert_((  "x0"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 1 )
    self.assert_( s.getBoundaryCondition() == False )
    s = m.getSpecies(3)
    self.assert_((  "x1"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == False )
    s = m.getSpecies(4)
    self.assert_((  "x2"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 1 )
    self.assert_( s.getBoundaryCondition() == False )
    s = m.getSpecies(5)
    self.assert_((  "x3"    == s.getName() ))
    self.assert_((  "cell"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( s.getBoundaryCondition() == False )
    self.assert_( m.getNumParameters() == 7 )
    p = m.getParameter(0)
    self.assert_((  "k1" == p.getName() ))
    self.assert_( p.getValue() == 1.2 )
    p = m.getParameter(1)
    self.assert_((  "k2" == p.getName() ))
    self.assert_( p.getValue() == 1000 )
    p = m.getParameter(2)
    self.assert_((  "k3" == p.getName() ))
    self.assert_( p.getValue() == 3000 )
    p = m.getParameter(3)
    self.assert_((  "k4" == p.getName() ))
    self.assert_( p.getValue() == 4.5 )
    self.assert_( m.getNumRules() == 4 )
    pr = m.getRule(0)
    self.assert_((  "t" == pr.getVariable() ))
    self.assert_((  "s1 + s2" == pr.getFormula() ))
    ud = pr.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 2 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_MOLE )
    self.assert_( ud.getUnit(0).getExponent() == 1 )
    self.assert_( ud.getUnit(1).getKind() == libsbml.UNIT_KIND_LITRE )
    self.assert_( ud.getUnit(1).getExponent() == -1 )
    self.assert_( pr.containsUndeclaredUnits() == False )
    pr = m.getRule(1)
    self.assert_((  "k" == pr.getVariable() ))
    self.assert_((  "k3/k2" == pr.getFormula() ))
    ud = pr.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 0 )
    self.assert_( pr.containsUndeclaredUnits() == True )
    scr = m.getRule(2)
    self.assert_((  "x2" == scr.getVariable() ))
    self.assert_((  "k * (s1+s2)/(1 + k)" == scr.getFormula() ))
    scr = m.getRule(3)
    self.assert_((  "x3" == scr.getVariable() ))
    self.assert_((  "p*(t - s2)" == scr.getFormula() ))
    self.assert_( m.getNumReactions() == 2 )
    r = m.getReaction(0)
    self.assert_((  "j1" == r.getName() ))
    self.assert_( r.getReversible() != False )
    self.assert_( r.getFast() == False )
    r = m.getReaction(1)
    self.assert_((  "j3" == r.getName() ))
    self.assert_( r.getReversible() != False )
    self.assert_( r.getFast() == False )
    r = m.getReaction(0)
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "x0" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "s1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    kl = r.getKineticLaw()
    self.assert_((  "k1 * x0" == kl.getFormula() ))
    r = m.getReaction(1)
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "s2" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "x1" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    kl = r.getKineticLaw()
    self.assert_((  "k4 * s2" == kl.getFormula() ))
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestReadFromFile3))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
