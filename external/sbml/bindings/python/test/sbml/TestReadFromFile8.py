#
# @file    TestReadFromFile8.py
# @brief   Reads test-data/l2v4-new.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Sarah Keating 
#
# $Id: TestReadFromFile8.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestReadFromFile8.py $
#
# This test file was converted from src/sbml/test/TestReadFromFile8.cpp
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

class TestReadFromFile8(unittest.TestCase):


  def test_read_l2v4_new(self):
    reader = libsbml.SBMLReader()
    filename = "../../sbml/test/test-data/"
    filename += "l2v4-new.xml"
    d = reader.readSBML(filename)
    if (d == None):
      pass    
    self.assert_( d.getLevel() == 2 )
    self.assert_( d.getVersion() == 4 )
    m = d.getModel()
    self.assert_( m != None )
    self.assert_( m.getId() ==  "l2v4_all" )
    self.assert_( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    self.assert_( c != None )
    self.assert_( c.getId() ==  "a" )
    self.assert_( c.getSize() == 1 )
    self.assertEqual( False, c.getConstant() )
    self.assert_( m.getNumEvents() == 1 )
    e = m.getEvent(0)
    self.assert_( e != None )
    self.assertEqual( True, e.getUseValuesFromTriggerTime() )
    self.assertEqual( True, e.isSetTrigger() )
    trigger = e.getTrigger()
    self.assert_( trigger != None )
    ast = trigger.getMath()
    self.assert_((  "lt(x, 3)" == libsbml.formulaToString(ast) ))
    self.assert_( e.getNumEventAssignments() == 1 )
    ea = e.getEventAssignment(0)
    self.assert_( ea != None )
    self.assert_( ea.getVariable() ==  "a" )
    ast = ea.getMath()
    self.assert_((  "x * p3" == libsbml.formulaToString(ast) ))
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestReadFromFile8))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
