/*
 *
 * @file    TestSBMLDocument.java
 * @brief   SBMLDocument unit tests
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Ben Bornstein 
 *
 * $Id: TestSBMLDocument.java 10124 2009-08-28 12:04:51Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/java/test/org/sbml/libsbml/test/sbml/TestSBMLDocument.java $
 *
 * This test file was converted from src/sbml/test/TestSBMLDocument.c
 * with the help of conversion sciprt (ctest_converter.pl).
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2009 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *--------------------------------------------------------------------------->*/


package org.sbml.libsbml.test.sbml;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestSBMLDocument {

  static void assertTrue(boolean condition) throws AssertionError
  {
    if (condition == true)
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      return;
    }
    else if ( (a == null) || (b == null) )
    {
      throw new AssertionError();
    }
    else if (a.equals(b))
    {
      return;
    }

    throw new AssertionError();
  }

  static void assertNotEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      throw new AssertionError();
    }
    else if ( (a == null) || (b == null) )
    {
      return;
    }
    else if (a.equals(b))
    {
      throw new AssertionError();
    }
  }

  static void assertEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(int a, int b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(int a, int b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }

  public void test_SBMLDocument_create()
  {
    SBMLDocument d = new  SBMLDocument();
    assertTrue( d.getTypeCode() == libsbml.SBML_DOCUMENT );
    assertTrue( d.getNotes() == null );
    assertTrue( d.getAnnotation() == null );
    assertTrue( d.getLevel() == 2 );
    assertTrue( d.getVersion() == 4 );
    assertTrue( d.getNumErrors() == 0 );
    d = null;
  }

  public void test_SBMLDocument_createWith()
  {
    SBMLDocument d = new  SBMLDocument(1,2);
    assertTrue( d.getTypeCode() == libsbml.SBML_DOCUMENT );
    assertTrue( d.getNotes() == null );
    assertTrue( d.getAnnotation() == null );
    assertTrue( d.getLevel() == 1 );
    assertTrue( d.getVersion() == 2 );
    assertTrue( d.getNumErrors() == 0 );
    d = null;
  }

  public void test_SBMLDocument_free_NULL()
  {
  }

  public void test_SBMLDocument_setLevelAndVersion()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,2,false);
    Model m1 = new  Model(2,2);
    d.setModel(m1);
    assertTrue( d.setLevelAndVersion(2,3,false) == true );
    assertTrue( d.setLevelAndVersion(2,1,false) == true );
    assertTrue( d.setLevelAndVersion(1,2,false) == true );
    assertTrue( d.setLevelAndVersion(1,1,false) == false );
    d = null;
  }

  public void test_SBMLDocument_setLevelAndVersion_Error()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,1,false);
    Model m1 = new  Model(2,1);
    Unit u = new  Unit(2,1);
    u.setKind(libsbml.UnitKind_forName("mole"));
    u.setOffset(3.2);
    UnitDefinition ud = new  UnitDefinition(2,1);
    ud.setId( "ud");
    ud.addUnit(u);
    m1.addUnitDefinition(ud);
    d.setModel(m1);
    assertTrue( d.setLevelAndVersion(2,2,false) == false );
    assertTrue( d.setLevelAndVersion(2,3,false) == false );
    assertTrue( d.setLevelAndVersion(1,2,false) == false );
    assertTrue( d.setLevelAndVersion(1,1,false) == false );
    d = null;
  }

  public void test_SBMLDocument_setLevelAndVersion_UnitsError()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,4,false);
    Model m1 = d.createModel();
    Compartment c = m1.createCompartment();
    c.setId( "c");
    Parameter p = m1.createParameter();
    p.setId( "p");
    p.setUnits( "mole");
    Rule r = m1.createAssignmentRule();
    r.setVariable( "c");
    r.setFormula( "p*p");
    assertTrue( d.setLevelAndVersion(2,2,false) == true );
    assertTrue( d.setLevelAndVersion(2,3,false) == true );
    assertTrue( d.setLevelAndVersion(1,2,false) == true );
    assertTrue( d.setLevelAndVersion(1,1,false) == false );
    d = null;
  }

  public void test_SBMLDocument_setLevelAndVersion_Warning()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,2,false);
    Model m1 = new  Model(2,2);
    (m1).setSBOTerm(2);
    d.setModel(m1);
    assertTrue( d.setLevelAndVersion(2,3,false) == true );
    assertTrue( d.setLevelAndVersion(2,1,false) == true );
    assertTrue( d.setLevelAndVersion(1,2,false) == true );
    assertTrue( d.setLevelAndVersion(1,1,false) == false );
    d = null;
  }

  public void test_SBMLDocument_setModel()
  {
    SBMLDocument d = new  SBMLDocument();
    Model m1 = new  Model(2,4);
    Model m2 = new  Model(2,4);
    assertTrue( d.getModel() == null );
    d.setModel(m1);
    assertTrue( !d.getModel().equals(m1) );
    d.setModel(d.getModel());
    assertTrue( !d.getModel().equals(m1) );
    d.setModel(m2);
    assertTrue( !d.getModel().equals(m2) );
    d = null;
  }

  public void test_SBMLDocument_setModel1()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,2,false);
    Model m1 = new  Model(2,1);
    long i = d.setModel(m1);
    assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
    assertTrue( d.getModel() == null );
    d = null;
  }

  public void test_SBMLDocument_setModel2()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,2,false);
    Model m1 = new  Model(1,2);
    long i = d.setModel(m1);
    assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
    assertTrue( d.getModel() == null );
    d = null;
  }

  public void test_SBMLDocument_setModel3()
  {
    SBMLDocument d = new  SBMLDocument();
    d.setLevelAndVersion(2,2,false);
    Model m1 = new  Model(2,2);
    long i = d.setModel(m1);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( d.getModel() != null );
    d = null;
  }

  /**
   * Loads the SWIG-generated libSBML Java module when this class is
   * loaded, or reports a sensible diagnostic message about why it failed.
   */
  static
  {
    String varname;
    String shlibname;

    if (System.getProperty("mrj.version") != null)
    {
      varname = "DYLD_LIBRARY_PATH";    // We're on a Mac.
      shlibname = "libsbmlj.jnilib and/or libsbml.dylib";
    }
    else
    {
      varname = "LD_LIBRARY_PATH";      // We're not on a Mac.
      shlibname = "libsbmlj.so and/or libsbml.so";
    }

    try
    {
      System.loadLibrary("sbmlj");
      // For extra safety, check that the jar file is in the classpath.
      Class.forName("org.sbml.libsbml.libsbml");
    }
    catch (SecurityException e)
    {
      e.printStackTrace();
      System.err.println("Could not load the libSBML library files due to a"+
                         " security exception.\n");
      System.exit(1);
    }
    catch (UnsatisfiedLinkError e)
    {
      e.printStackTrace();
      System.err.println("Error: could not link with the libSBML library files."+
                         " It is likely\nyour " + varname +
                         " environment variable does not include the directories\n"+
                         "containing the " + shlibname + " library files.\n");
      System.exit(1);
    }
    catch (ClassNotFoundException e)
    {
      e.printStackTrace();
      System.err.println("Error: unable to load the file libsbmlj.jar."+
                         " It is likely\nyour -classpath option and CLASSPATH" +
                         " environment variable\n"+
                         "do not include the path to libsbmlj.jar.\n");
      System.exit(1);
    }
  }
}
