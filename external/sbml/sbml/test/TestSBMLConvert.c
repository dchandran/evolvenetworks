/**
 * \file    TestSBMLConvert.c
 * \brief   SBMLConvert unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestSBMLConvert.c 10129 2009-08-28 12:23:22Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestSBMLConvert.c $
 */
/* Copyright 2003 California Institute of Technology and
 * Japan Science and Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * documentation provided hereunder is on an "as is" basis, and the
 * California Institute of Technology and Japan Science and Technology
 * Corporation have no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall the
 * California Institute of Technology or the Japan Science and Technology
 * Corporation be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if the
 * California Institute of Technology and/or Japan Science and Technology
 * Corporation have been advised of the possibility of such damage.  See
 * the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * The original code contained here was initially developed by:
 *
 *     Ben Bornstein
 *     The Systems Biology Markup Language Development Group
 *     ERATO Kitano Symbiotic Systems Project
 *     Control and Dynamical Systems, MC 107-81
 *     California Institute of Technology
 *     Pasadena, CA, 91125, USA
 *
 *     http://www.cds.caltech.edu/erato
 *     mailto:sbml-team@caltech.edu
 *
 * Contributor(s):
 */


#include <sbml/common/common.h>
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLTypes.h>

#include <check.h>


START_TEST (test_SBMLConvert_convertToL2_SBMLDocument)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(1, 2);


  SBMLDocument_setLevelAndVersion(d, 2, 1);

  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 1, NULL );

  SBMLDocument_setLevelAndVersion(d, 2, 2);
  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 2, NULL );

  SBMLDocument_setLevelAndVersion(d, 2, 3);
  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 3, NULL );

  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_addModifiersToReaction)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(1, 2);
  Model_t        *m = SBMLDocument_createModel(d);
  Reaction_t     *r = Model_createReaction(m);
  KineticLaw_t  *kl = Reaction_createKineticLaw(r);
  KineticLaw_setFormula(kl, "k1*S1*S2*S3*S4*S5");

  SpeciesReference_t *ssr1;
  SpeciesReference_t *ssr2;


  Species_t *s1 = Model_createSpecies( m ); 
  Species_setId( s1, "S1" );
  Species_t *s2 = Model_createSpecies( m ); 
  Species_setId( s2, "S2");
  Species_t *s3 = Model_createSpecies( m ); 
  Species_setId( s3, "S3");
  Species_t *s4 = Model_createSpecies( m ); 
  Species_setId( s4, "S4");
  Species_t *s5 = Model_createSpecies( m ); 
  Species_setId( s5, "S5");

  SpeciesReference_t *sr1 = Reaction_createReactant( r );
  SpeciesReference_t *sr2 = Reaction_createReactant( r );
  SpeciesReference_t *sr3 = Reaction_createProduct ( r );

  SpeciesReference_setSpecies(sr1, "S1");
  SpeciesReference_setSpecies(sr2, "S2");
  SpeciesReference_setSpecies(sr3, "S5");

  fail_unless( Reaction_getNumModifiers(r) == 0, NULL );

  SBMLDocument_setLevelAndVersion(d, 2, 1);

  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 1, NULL );


  fail_unless( Reaction_getNumModifiers(Model_getReaction(m, 0)) == 2, NULL );

  ssr1 = (SpeciesReference_t *) Reaction_getModifier(Model_getReaction(m, 0), 0);
  ssr2 = (SpeciesReference_t *) Reaction_getModifier(Model_getReaction(m, 0), 1);

  fail_unless( !strcmp(SpeciesReference_getSpecies(ssr1), "S3"), NULL );
  fail_unless( !strcmp(SpeciesReference_getSpecies(ssr2), "S4"), NULL );

  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_convertToL1_SBMLDocument)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(2, 1);

  SBMLDocument_setLevelAndVersion(d, 1, 2);

  fail_unless( SBMLDocument_getLevel  (d) == 1, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 2, NULL );

  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_convertToL1_Species_Amount)
{
  SBMLDocument_t *d   = SBMLDocument_createWithLevelAndVersion(2, 1);
  Model_t        *m   = SBMLDocument_createModel(d);
  const char     *sid = "C";
  Compartment_t  *c   = Compartment_create(2, 4);
  Species_t      *s   = Species_create(2, 4);


  Compartment_setId   ( c, sid );
  Model_addCompartment( m, c   );

  Species_setCompartment  ( s, sid  ); 
  Species_setInitialAmount( s, 2.34 );
  Model_addSpecies        ( m, s    );
  
  SBMLDocument_setLevelAndVersion(d, 1, 2);

  fail_unless( Species_getInitialAmount(s) == 2.34, NULL );

  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_convertToL1_Species_Concentration)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(2, 1);
  Model_t        *m = SBMLDocument_createModel(d);
  const char   *sid = "C";
  Compartment_t  *c = 
    Compartment_create(2, 1);
  Species_t      *s = 
    Species_create(2, 1);


  Compartment_setId   ( c, sid );
  Compartment_setSize ( c, 1.2 ); 
  Model_addCompartment( m, c   );

  Species_setId                  ( s, "s"  );
  Species_setCompartment         ( s, sid  ); 
  Species_setInitialConcentration( s, 2.34 );
  Model_addSpecies               ( m, s    );
  
  SBMLDocument_setLevelAndVersion(d, 1, 2);

  /**
   * These tests will fail under Cygwin because of a minimal
   * setlocale() implementation (see setlocale manpage).
   */
#ifndef CYGWIN
  fail_unless( Species_getInitialAmount(Model_getSpecies(m, 0)) == 2.808, NULL );
#endif

  Species_t * s1 = Model_getSpecies(m, 0);
  fail_unless (s1 != NULL);
  fail_unless (!strcmp(Species_getCompartment(s1), "C"));
  fail_unless(Compartment_getSize(Model_getCompartmentById(m, "C")) == 1.2);
  fail_unless(Species_getInitialConcentration(s1) == 2.34);
  fail_unless(Species_isSetInitialConcentration(s1) == 1);

  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_convertToL2v4_DuplicateAnnotations_doc)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(2, 1);
  SBMLDocument_createModel(d);

  char * annotation = "<rdf/>\n<rdf/>";

  int i = SBase_setAnnotationString((SBase_t *) (d), annotation);
  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 1, NULL );
  fail_unless( XMLNode_getNumChildren(SBase_getAnnotation((SBase_t *) (d))) == 2);

  SBMLDocument_setLevelAndVersion(d, 2, 4);

  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 4, NULL );

  fail_unless( XMLNode_getNumChildren(SBase_getAnnotation((SBase_t *) (d))) == 1);


  SBMLDocument_free(d);
}
END_TEST


START_TEST (test_SBMLConvert_convertToL2v4_DuplicateAnnotations_model)
{
  SBMLDocument_t *d = SBMLDocument_createWithLevelAndVersion(2, 1);
  Model_t * m = SBMLDocument_createModel(d);

  char * annotation = "<rdf/>\n<rdf/>";

  int i = SBase_setAnnotationString((SBase_t *) (m), annotation);
  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 1, NULL );
  fail_unless( XMLNode_getNumChildren(SBase_getAnnotation((SBase_t *) (m))) == 2);

  SBMLDocument_setLevelAndVersion(d, 2, 4);

  fail_unless( SBMLDocument_getLevel  (d) == 2, NULL );
  fail_unless( SBMLDocument_getVersion(d) == 4, NULL );

  m = SBMLDocument_getModel(d);
  fail_unless( XMLNode_getNumChildren(SBase_getAnnotation((SBase_t *) (m))) == 1);


  SBMLDocument_free(d);
}
END_TEST


Suite *
create_suite_SBMLConvert (void) 
{ 
  Suite *suite = suite_create("SBMLConvert");
  TCase *tcase = tcase_create("SBMLConvert");


  tcase_add_test( tcase, test_SBMLConvert_convertToL2_SBMLDocument       );
  tcase_add_test( tcase, test_SBMLConvert_addModifiersToReaction         );
  tcase_add_test( tcase, test_SBMLConvert_convertToL1_SBMLDocument          );
  tcase_add_test( tcase, test_SBMLConvert_convertToL1_Species_Amount        );
  tcase_add_test( tcase, test_SBMLConvert_convertToL1_Species_Concentration );
  tcase_add_test( tcase, test_SBMLConvert_convertToL2v4_DuplicateAnnotations_doc );
  tcase_add_test( tcase, test_SBMLConvert_convertToL2v4_DuplicateAnnotations_model );

  suite_add_tcase(suite, tcase);

  return suite;
}
