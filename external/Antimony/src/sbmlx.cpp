#include <cassert>
#include <string>
#include <vector>

extern bool CaselessStrCmp(const std::string& lhs, const std::string& rhs);

#ifndef NSBML
#include "sbmlx.h"
#include "formula.h"
#include "variable.h"
#include "registry.h"

using namespace std;

// SBase objects no longer have IDs :( //update: now they do again!
string getNameFromSBMLObject(const SBase* sbml, string basename)
{
  string name = sbml->getId();
  if (name == "") {
    name = sbml->getName();
    //Names can have spaces, so...
    while (name.find(" ") != string::npos) {
      name.replace(name.find(" "), 1, "_");
    }
  }
  if (name=="") {
    long num=0;
    Variable* foundvar = NULL;
    do {
      char charnum[50];
      sprintf(charnum, "%li", num);
      num++;
      name = basename;
      name += charnum;
      vector<string> fullname;
      fullname.push_back(name);
      foundvar = g_registry.CurrentModule()->GetVariable(fullname);
    } while (foundvar != NULL);
  }
  assert(name != "");
  return name;
}

/*
string getNameFromSBMLObject(string ID, string name, string basename)
{
  if (ID != "") return ID;
  if (name != "") return name;
  long num=0;
  Variable* foundvar = NULL;
  do {
    char charnum[50];
    sprintf(charnum, "%li", num);
    num++;
    name = basename;
    name += charnum;
    vector<string> fullname;
    fullname.push_back(name);
    foundvar = g_registry.CurrentModule()->GetVariable(fullname);
  } while (foundvar != NULL);
  assert(name != "");
  return name;
}
*/
void setTimeName(ASTNode *node)  
{
  if (node->getType() == AST_NAME_TIME) {
    node->setName("time");
  }
  for (unsigned int c = 0; c < node->getNumChildren() ; c++) {
    setTimeName(node->getChild(c));
  }
}

string parseASTNodeToString(const ASTNode* ASTform) {
  if (ASTform==NULL) return "";
  ASTNode clone(*ASTform);
  setTimeName(&clone);
  char* formula = SBML_formulaToString(&clone);
  string ret = formula;
#ifndef WIN32
  free(formula);
#endif
  return ret;
}


void setTimeType(ASTNode_t* node)  
{
  if (node->isOperator() == false && node->isNumber() == false) {
    if (string(node->getName()) == "time") {
      node->setType(AST_NAME_TIME);
    }
  }
  for (unsigned int c = 0; c < node->getNumChildren() ; c++) {
    setTimeType(node->getChild(c));
  }
}

ASTNode* parseStringToASTNode(const string& formula) {
  ASTNode* rootnode = SBML_parseFormula(formula.c_str());
  if (rootnode == NULL) return NULL;
  if (formula.find("time") != string::npos) {
    setTimeType(rootnode);
  }
  return rootnode;
}

#endif


//SBML models might have variable names in them that are reserved keywords in Antimony (like 'compartment', to take a huge example).  FixName fixes this so that you can output readable Antimony again.
bool FixName(std::string& name)
{
  while (name.size() && name[0] == ' ') {
    name.erase(0, 1);
  }
  while (name.size() && name[name.size()-1] == ' ') {
    name.erase(name.size()-1, 1);
  }
  const char* keywords[] = {
  "DNA",
  "at",
  "compartment",
  "const",
  "end",
  "event",
  "ext",
  "formula",
  "function",
  "gene",
  "import",
  "in",
  "is",
  "model",
  "module",
  "operator",
  "reaction",
  "species",
  "var",

  "nan", //not a number in Math

  "abs"
  , "and"
  , "annotation"
  , "annotation-xml"
  , "apply"
  , "arccos"
  , "arccosh"
  , "arccot"
  , "arccoth"
  , "arccsc"
  , "arccsch"
  , "arcsec"
  , "arcsech"
  , "arcsin"
  , "arcsinh"
  , "arctan"
  , "arctanh"
  , "bvar"
  , "ceiling"
  //  , "ci"
  , "cn"
  , "cos"
  , "cosh"
  , "cot"
  , "coth"
  , "csc"
  , "csch"
  , "csymbol"
  , "degree"
  , "divide"
  , "eq"
  , "exp"
  , "exponentiale"
  , "factorial"
  , "false"
  , "floor"
  , "geq"
  , "gt"
  , "infinity"
  , "lambda"
  , "leq"
  , "ln"
  , "log"
  , "logbase"
  , "lt"
  , "math"
  , "minus"
  , "neq"
  , "not"
  , "notanumber"
  , "or"
  , "otherwise"
  , "pi"
  , "piece"
  , "piecewise"
  , "plus"
  , "power"
  , "root"
  , "sec"
  , "sech"
  , "semantics"
  , "sep"
  , "sin"
  , "sinh"
  , "tan"
  , "tanh"
  , "times"
  , "true"
  , "xor"
  , "acos"
  , "asin"
  , "atan"
  , "ceil"
  , "delay"
  , "log10"
  , "pow"
  , "sqr"
  , "sqrt"
  , "time"
  };
  for (size_t kw=0; kw<98; kw++) {
    if (CaselessStrCmp(name, keywords[kw])) {
      name += "_";
      return true;
    }
  }
  for (size_t pos=0; pos<name.size(); pos++) {
    if (!(isalpha(name[pos]) || isdigit(name[pos]) || name[pos]=='_')) {
      name[pos] = '_';
    }
  }
  return false;
}

bool FixName(vector<string>& names)
{
  bool ret = false;
  for (size_t n=0; n<names.size(); n++) {
    if (FixName(names[n])) {
      ret = true;
    }
  }
  return ret;
}

void FixName(vector<vector<string> >& allnames) 
{
  for (size_t n=0; n<allnames.size(); n++) {
    FixName(allnames[n]);
  }
}


void FixName(map<vector<string>, Variable*>& varmap)
{
  for (map<vector<string>, Variable*>::iterator vm=varmap.begin(); vm != varmap.end();)
  {
    vector<string> name = vm->first;
    if (FixName(name)) {
      map<vector<string>, Variable*>::iterator vm2 = vm;
      vm++;
      varmap.insert(make_pair(name, vm2->second));
      varmap.erase(vm2);
    }
    else {
      vm++;
    }
  }
  
}
