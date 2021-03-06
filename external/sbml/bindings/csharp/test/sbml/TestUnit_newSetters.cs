/// 
///  @file    TestUnit_newSetters.cs
///  @brief   Unit unit tests for new set function API
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id$
///  $HeadURL$
/// 
///  This test file was converted from src/sbml/test/TestUnit_newSetters.c
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

  public class TestUnit_newSetters {
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
      else if ( (a == null) || (b == null) )
      {
        throw new AssertionError();
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
      else if ( (a == null) || (b == null) )
      {
        return;
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

    private Unit U;

    public void setUp()
    {
      U = new  Unit(1,2);
      if (U == null);
      {
      }
    }

    public void tearDown()
    {
      U = null;
    }

    public void test_Unit_removeScale()
    {
      long i = U.setScale(2);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( U.getScale() == 2 );
      i = Unit.removeScale(U);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( U.getScale() == 0 );
      assertTrue( U.getMultiplier() == 100 );
    }

    public void test_Unit_setExponent1()
    {
      long i = U.setExponent(2);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( U.getExponent() == 2 );
    }

    public void test_Unit_setKind1()
    {
      long i = U.setKind(libsbml.UnitKind_forName("cell"));
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      assertEquals( false, U.isSetKind() );
    }

    public void test_Unit_setKind2()
    {
      long i = U.setKind(libsbml.UnitKind_forName("litre"));
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( true, U.isSetKind() );
    }

    public void test_Unit_setMultiplier1()
    {
      long i = U.setMultiplier(2);
      assertTrue( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE );
      assertTrue( U.getMultiplier() == 2 );
    }

    public void test_Unit_setMultiplier2()
    {
      Unit c = new  Unit(2,2);
      long i = c.setMultiplier(4);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( c.getMultiplier() == 4 );
      c = null;
    }

    public void test_Unit_setOffset1()
    {
      long i = U.setOffset(2.0);
      assertTrue( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE );
      assertTrue( U.getOffset() == 0 );
    }

    public void test_Unit_setOffset2()
    {
      Unit U1 = new  Unit(2,1);
      long i = U1.setOffset(2.0);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( U1.getOffset() == 2 );
    }

    public void test_Unit_setScale1()
    {
      long i = U.setScale(2);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( U.getScale() == 2 );
    }

  }
}
