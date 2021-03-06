#
# @file    TestRateRule.py
# @brief   RateRule unit tests
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestRateRule.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestRateRule.py $
#
# This test file was converted from src/sbml/test/TestRateRule.c
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

class TestRateRule(unittest.TestCase):

  RR = None

  def setUp(self):
    self.RR = libsbml.RateRule(1,2)
    if (self.RR == None):
      pass    
    pass  

  def tearDown(self):
    self.RR = None
    pass  

  def test_RateRule_create(self):
    self.assert_( self.RR.getTypeCode() == libsbml.SBML_RATE_RULE )
    self.assert_( self.RR.getMetaId() == "" )
    self.assert_( self.RR.getNotes() == None )
    self.assert_( self.RR.getAnnotation() == None )
    self.assert_( self.RR.getFormula() == "" )
    self.assert_( self.RR.getMath() == None )
    self.assert_( self.RR.getVariable() == "" )
    self.assert_( self.RR.getType() == libsbml.RULE_TYPE_RATE )
    pass  

  def test_RateRule_createWithNS(self):
    xmlns = libsbml.XMLNamespaces()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = libsbml.SBMLNamespaces(2,1)
    sbmlns.addNamespaces(xmlns)
    object = libsbml.RateRule(sbmlns)
    self.assert_( object.getTypeCode() == libsbml.SBML_RATE_RULE )
    self.assert_( object.getMetaId() == "" )
    self.assert_( object.getNotes() == None )
    self.assert_( object.getAnnotation() == None )
    self.assert_( object.getLevel() == 2 )
    self.assert_( object.getVersion() == 1 )
    self.assert_( object.getNamespaces() != None )
    self.assert_( object.getNamespaces().getLength() == 2 )
    object = None
    pass  

  def test_RateRule_free_NULL(self):
    pass  

  def test_RateRule_setVariable(self):
    variable =  "x";
    self.RR.setVariable(variable)
    self.assert_(( variable == self.RR.getVariable() ))
    self.assertEqual( True, self.RR.isSetVariable() )
    if (self.RR.getVariable() == variable):
      pass    
    self.RR.setVariable(self.RR.getVariable())
    self.assert_(( variable == self.RR.getVariable() ))
    self.RR.setVariable("")
    self.assertEqual( False, self.RR.isSetVariable() )
    if (self.RR.getVariable() != None):
      pass    
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestRateRule))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
