#
# @file    TestReadFromFile4.py
# @brief   Reads tests/l1v1-minimal.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestReadFromFile4.py 8704 2009-01-04 02:26:05Z mhucka $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestReadFromFile4.py $
#
# This test file was converted from src/sbml/test/TestReadFromFile4.c
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

class TestReadFromFile4(unittest.TestCase):


  def test_read_l1v1_minimal(self):
    filename = "../../sbml/test/test-data/l1v1-minimal.xml"
    d = libsbml.readSBML(filename)
    if (d == None):
      pass    
    self.assert_( d.getLevel() == 1 )
    self.assert_( d.getVersion() == 1 )
    m = d.getModel()
    self.assert_( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    self.assert_((  "x" == c.getName() ))
    self.assert_( m.getNumSpecies() == 1 )
    s = m.getSpecies(0)
    self.assert_((  "y"  == s.getName() ))
    self.assert_((  "x"  == s.getCompartment() ))
    self.assert_( s.getInitialAmount() == 1 )
    self.assert_( s.getBoundaryCondition() == False )
    self.assert_( m.getNumReactions() == 1 )
    r = m.getReaction(0)
    self.assert_((  "r" == r.getName() ))
    self.assert_( r.getReversible() != False )
    self.assert_( r.getFast() == False )
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 1 )
    sr = r.getReactant(0)
    self.assert_((  "y" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    sr = r.getProduct(0)
    self.assert_((  "y" == sr.getSpecies() ))
    self.assert_( sr.getStoichiometry() == 1 )
    self.assert_( sr.getDenominator() == 1 )
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestReadFromFile4))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
