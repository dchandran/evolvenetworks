/*
 *
 * @file    TestValidASTNode.java
 * @brief   Test the isWellFormedASTNode function
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Sarah Keating
 *
 * $Id$
 * $HeadURL$
 *
 * This test file was converted from src/sbml/test/TestValidASTNode.cpp
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


package org.sbml.libsbml.test.math;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestValidASTNode {

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

  public void test_ValidASTNode_Name()
  {
    ASTNode n = libsbml.parseFormula("c");
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( false, (n.isWellFormedASTNode()) );
    n = null;
  }

  public void test_ValidASTNode_Number()
  {
    ASTNode n = libsbml.parseFormula("1.2");
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( false, (n.isWellFormedASTNode()) );
    n = null;
  }

  public void test_ValidASTNode_binary()
  {
    ASTNode n = new ASTNode(libsbml.AST_DIVIDE);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode c = libsbml.parseFormula("c");
    n.addChild(c);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( true, n.isWellFormedASTNode() );
    n = null;
  }

  public void test_ValidASTNode_lambda()
  {
    ASTNode n = new ASTNode(libsbml.AST_LAMBDA);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode c = libsbml.parseFormula("c");
    n.addChild(c);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode e = libsbml.parseFormula("e");
    n.addChild(e);
    assertEquals( true, n.isWellFormedASTNode() );
    n = null;
  }

  public void test_ValidASTNode_nary()
  {
    ASTNode n = new ASTNode(libsbml.AST_TIMES);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode c = libsbml.parseFormula("c");
    n.addChild(c);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode e = libsbml.parseFormula("e");
    n.addChild(e);
    assertEquals( true, n.isWellFormedASTNode() );
    n = null;
  }

  public void test_ValidASTNode_root()
  {
    ASTNode n = new ASTNode(libsbml.AST_FUNCTION_ROOT);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode c = libsbml.parseFormula("c");
    n.addChild(c);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode d = libsbml.parseFormula("3");
    n.addChild(d);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode e = libsbml.parseFormula("3");
    n.addChild(e);
    assertEquals( false, (n.isWellFormedASTNode()) );
    n = null;
  }

  public void test_ValidASTNode_setType()
  {
    ASTNode n = new ASTNode();
    long i = n.setType(libsbml.AST_REAL);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( n.getType() == libsbml.AST_REAL );
    i = n.setType(libsbml.AST_PLUS);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( n.getType() == libsbml.AST_PLUS );
    char c = n.getCharacter();
    assertTrue(c == '+');
    i = n.setType(libsbml.AST_FUNCTION_ARCCOSH);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( n.getType() == libsbml.AST_FUNCTION_ARCCOSH );
    i = n.setType(libsbml.AST_UNKNOWN);
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    assertTrue( n.getType() == libsbml.AST_UNKNOWN );
    n = null;
  }

  public void test_ValidASTNode_unary()
  {
    ASTNode n = new ASTNode(libsbml.AST_FUNCTION_ABS);
    assertEquals( false, (n.isWellFormedASTNode()) );
    ASTNode c = libsbml.parseFormula("c");
    n.addChild(c);
    assertEquals( true, n.isWellFormedASTNode() );
    ASTNode d = libsbml.parseFormula("d");
    n.addChild(d);
    assertEquals( false, (n.isWellFormedASTNode()) );
    n = null;
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