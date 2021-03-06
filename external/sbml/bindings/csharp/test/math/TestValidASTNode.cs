/// 
///  @file    TestValidASTNode.cs
///  @brief   Test the isWellFormedASTNode function
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id$
///  $HeadURL$
/// 
///  This test file was converted from src/sbml/test/TestValidASTNode.cpp
///  with the help of conversion sciprt (ctest_converter.pl).
/// 
/// <!---------------------------------------------------------------------------
///  This file is part of libSBML.  Please visit http://sbml.org for more
///  information about SBML, and the latest version of libSBML.
/// 
///  Copyright 2005-2009 California Institute of Technology.
///  Copyright 2002-2005 California Institute of Technology and
///                      Japan Science and Technology Corporation.
///  
///  This library is free software; you can redistribute it and/or modify it
///  under the terms of the GNU Lesser General Public License as published by
///  the Free Software Foundation.  A copy of the license agreement is provided
///  in the file named "LICENSE.txt" included with this software distribution
///  and also available online as http://sbml.org/software/libsbml/license.html
/// --------------------------------------------------------------------------->*/


namespace LibSBMLCSTest {

  using libsbml;

  using  System.IO;

  public class TestValidASTNode {
    public class AssertionError : System.Exception 
    {
      public AssertionError() : base()
      {
        
      }
    }


    static void assertTrue(bool condition)
    {
      if (condition == true)
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        return;
      }
      else if (a.Equals(b))
      {
        return;
      }
  
      throw new AssertionError();
    }

    static void assertNotEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        throw new AssertionError();
      }
      else if (a.Equals(b))
      {
        throw new AssertionError();
      }
    }

    static void assertEquals(bool a, bool b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(bool a, bool b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(int a, int b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(int a, int b)
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
      assertTrue(('+'  == c ));
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

  }
}
