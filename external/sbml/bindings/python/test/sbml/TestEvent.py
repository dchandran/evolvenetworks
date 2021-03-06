#
# @file    TestEvent.py
# @brief   SBML Event unit tests
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestEvent.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestEvent.py $
#
# This test file was converted from src/sbml/test/TestEvent.c
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

class TestEvent(unittest.TestCase):

  E = None

  def setUp(self):
    self.E = libsbml.Event(2,4)
    if (self.E == None):
      pass    
    pass  

  def tearDown(self):
    self.E = None
    pass  

  def test_Event_create(self):
    self.assert_( self.E.getTypeCode() == libsbml.SBML_EVENT )
    self.assert_( self.E.getMetaId() == "" )
    self.assert_( self.E.getNotes() == None )
    self.assert_( self.E.getAnnotation() == None )
    self.assert_( self.E.getId() == "" )
    self.assert_( self.E.getName() == "" )
    self.assert_( self.E.getTrigger() == None )
    self.assert_( self.E.getDelay() == None )
    self.assert_( self.E.getTimeUnits() == "" )
    self.assert_( self.E.getNumEventAssignments() == 0 )
    pass  

  def test_Event_createWithNS(self):
    xmlns = libsbml.XMLNamespaces()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = libsbml.SBMLNamespaces(2,4)
    sbmlns.addNamespaces(xmlns)
    object = libsbml.Event(sbmlns)
    self.assert_( object.getTypeCode() == libsbml.SBML_EVENT )
    self.assert_( object.getMetaId() == "" )
    self.assert_( object.getNotes() == None )
    self.assert_( object.getAnnotation() == None )
    self.assert_( object.getLevel() == 2 )
    self.assert_( object.getVersion() == 4 )
    self.assert_( object.getNamespaces() != None )
    self.assert_( object.getNamespaces().getLength() == 2 )
    object = None
    pass  

  def test_Event_free_NULL(self):
    pass  

  def test_Event_full(self):
    math1 = libsbml.parseFormula("0")
    trigger = libsbml.Trigger(2,4)
    math = libsbml.parseFormula("0")
    e = libsbml.Event(2,4)
    ea = libsbml.EventAssignment(2,4)
    ea.setVariable( "k")
    ea.setMath(math)
    trigger.setMath(math1)
    e.setTrigger(trigger)
    e.setId( "e1")
    e.setName( "Set k2 to zero when P1 <= t")
    e.addEventAssignment(ea)
    self.assert_( e.getNumEventAssignments() == 1 )
    self.assert_( e.getEventAssignment(0) != ea )
    math = None
    e = None
    pass  

  def test_Event_removeEventAssignment(self):
    o1 = self.E.createEventAssignment()
    o2 = self.E.createEventAssignment()
    o3 = self.E.createEventAssignment()
    o3.setVariable("test")
    self.assert_( self.E.removeEventAssignment(0) == o1 )
    self.assert_( self.E.getNumEventAssignments() == 2 )
    self.assert_( self.E.removeEventAssignment(0) == o2 )
    self.assert_( self.E.getNumEventAssignments() == 1 )
    self.assert_( self.E.removeEventAssignment("test") == o3 )
    self.assert_( self.E.getNumEventAssignments() == 0 )
    o1 = None
    o2 = None
    o3 = None
    pass  

  def test_Event_setDelay(self):
    math1 = libsbml.parseFormula("0")
    Delay = libsbml.Delay(2,4)
    Delay.setMath(math1)
    self.E.setDelay(Delay)
    self.assert_( self.E.getDelay() != None )
    self.assertEqual( True, self.E.isSetDelay() )
    if (self.E.getDelay() == Delay):
      pass    
    self.E.setDelay(self.E.getDelay())
    self.assert_( self.E.getDelay() != Delay )
    self.E.setDelay(None)
    self.assertEqual( False, self.E.isSetDelay() )
    if (self.E.getDelay() != None):
      pass    
    pass  

  def test_Event_setId(self):
    id =  "e1";
    self.E.setId(id)
    self.assert_(( id == self.E.getId() ))
    self.assertEqual( True, self.E.isSetId() )
    if (self.E.getId() == id):
      pass    
    self.E.setId(self.E.getId())
    self.assert_(( id == self.E.getId() ))
    self.E.setId("")
    self.assertEqual( False, self.E.isSetId() )
    if (self.E.getId() != None):
      pass    
    pass  

  def test_Event_setName(self):
    name =  "Set_k2";
    self.E.setName(name)
    self.assert_(( name == self.E.getName() ))
    self.assertEqual( True, self.E.isSetName() )
    if (self.E.getName() == name):
      pass    
    self.E.setName(self.E.getName())
    self.assert_(( name == self.E.getName() ))
    self.E.setName("")
    self.assertEqual( False, self.E.isSetName() )
    if (self.E.getName() != None):
      pass    
    pass  

  def test_Event_setTimeUnits(self):
    E1 = libsbml.Event(2,1)
    units =  "second";
    E1.setTimeUnits(units)
    self.assert_(( units == E1.getTimeUnits() ))
    self.assertEqual( True, E1.isSetTimeUnits() )
    if (E1.getTimeUnits() == units):
      pass    
    E1.setTimeUnits(E1.getTimeUnits())
    self.assert_(( units == E1.getTimeUnits() ))
    E1.setTimeUnits("")
    self.assertEqual( False, E1.isSetTimeUnits() )
    if (E1.getTimeUnits() != None):
      pass    
    E1 = None
    pass  

  def test_Event_setTrigger(self):
    math1 = libsbml.parseFormula("0")
    trigger = libsbml.Trigger(2,4)
    trigger.setMath(math1)
    self.E.setTrigger(trigger)
    self.assert_( self.E.getTrigger() != None )
    self.assertEqual( True, self.E.isSetTrigger() )
    if (self.E.getTrigger() == trigger):
      pass    
    self.E.setTrigger(self.E.getTrigger())
    self.assert_( self.E.getTrigger() != trigger )
    self.E.setTrigger(None)
    self.assertEqual( False, self.E.isSetTrigger() )
    if (self.E.getTrigger() != None):
      pass    
    pass  

  def test_Event_setUseValuesFromTriggerTime(self):
    object = libsbml.Event(2,4)
    object.setUseValuesFromTriggerTime(False)
    self.assert_( object.getUseValuesFromTriggerTime() == False )
    object.setUseValuesFromTriggerTime(True)
    self.assert_( object.getUseValuesFromTriggerTime() == True )
    object = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestEvent))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
