#ifndef NCELLML

#include "cellmlx.h"
#include "MathMLInputServices.h"
#include "ICellMLInputServices.h"
#include "IAnnoTools.h"
#include <nsStringAPI.h>
#include <nsDebug.h>
#include <nsCOMPtr.h>
#include <prtypes.h>
#include <nsServiceManagerUtils.h>

#define CELLML_BOOTSTRAP_CONTRACTID "@cellml.org/cellml-bootstrap;1"

void Module::LoadCellMLModel(nsCOMPtr<cellml_apiIModel> model, vector<nsCOMPtr<cellml_apiICellMLComponent> > top_components)
{
  assert(m_cellmlcomponent==NULL);
  if (m_cellmlmodel != NULL) {
    assert(false); 
  }
  m_cellmlmodel = model;

  //Translate
  nsString cellmltext;

  //Components become sub-modules
  for (size_t comp=0; comp<top_components.size(); comp++) {
    //Load the component as a submodule (which should already have been loaded by the registry in LoadCellML)
    string compname = GetNameAccordingToEncapsulationParent(top_components[comp], model);
    string cellmlmodname = GetModuleNameFrom(top_components[comp]);
    Variable* var = AddOrFindVariable(&compname);
    if(var->SetModule(&cellmlmodname)) {
      assert(false);
      return;
    }
  }

  FixNames();  //In case the name of one of the modules is something like 'time'.
}

bool HasTimeUnits(nsCOMPtr<cellml_apiICellMLVariable> cmlvar) {
  nsresult rv;
  nsString cellmltext;
  rv = cmlvar->GetUnitsName(cellmltext);
  if (cellmltext == ToNSString("second")) return true;
  nsCOMPtr<cellml_apiIUnits> units;
  rv = cmlvar->GetUnitsElement(getter_AddRefs(units));
  if (units==NULL) return false;
  nsCOMPtr<cellml_apiIUnitSet> unitset;
  rv = units->GetUnitCollection(getter_AddRefs(unitset));
  nsCOMPtr<cellml_apiIUnitIterator> uniti;
  rv = unitset->IterateUnits(getter_AddRefs(uniti));
  nsCOMPtr<cellml_apiIUnit> unit;
  rv = uniti->NextUnit(getter_AddRefs(unit));
  while (unit != NULL) {
    rv = unit->GetUnits(cellmltext);
    if (cellmltext == ToNSString("second")) return true;
    rv = uniti->NextUnit(getter_AddRefs(unit));
  }
  return false;
}

void Module::LoadCellMLComponent(nsCOMPtr<cellml_apiICellMLComponent> component)
{
  nsString cellmltext;
  string cellmlname;
  nsresult rv;

  assert(m_cellmlmodel==NULL);
  if (m_cellmlcomponent != NULL) {
    assert(false);
  }
  //Variables
  nsCOMPtr<cellml_apiICellMLVariableSet> varset;
  rv = component->GetVariables(getter_AddRefs(varset));
  nsCOMPtr<cellml_apiICellMLVariableIterator> vsi;
  rv = varset->IterateVariables(getter_AddRefs(vsi));
  nsCOMPtr<cellml_apiICellMLVariable> cmlvar;
  rv = vsi->NextVariable(getter_AddRefs(cmlvar));
  while (cmlvar != NULL) {
    rv = cmlvar->GetName(cellmltext);
    cellmlname = ToThinString(cellmltext.get());
    Variable* antvar = AddOrFindVariable(&cellmlname);
    //antvar->SetIsConst(false); //This forces it to be output in the Antimony script even if it's not otherwise used.
    rv = cmlvar->GetInitialValue(cellmltext);
    cellmlname = ToThinString(cellmltext.get());
    if (cellmlname != "") {
      Formula* formula = g_registry.NewBlankFormula();
      setFormulaWithString(cellmlname, formula, this);
      antvar->SetFormula(formula);
    }
    //Put it in the module interface if it has one:
    PRUint32 vi;
    rv = cmlvar->GetPublicInterface(&vi);
    if (vi != 2) {
      AddVariableToExportList(antvar);
    }
    else {
      rv = cmlvar->GetPrivateInterface(&vi);
      if (vi != 2) {
        AddVariableToExportList(antvar);
      }
    }
    
    rv = vsi->NextVariable(getter_AddRefs(cmlvar));
  }

  //Reactions
  nsCOMPtr<cellml_apiIReactionSet> rxnset;
  rv = component->GetReactions(getter_AddRefs(rxnset));
  nsCOMPtr<cellml_apiIReactionIterator> rxni;
  rv = rxnset->IterateReactions(getter_AddRefs(rxni));
  nsCOMPtr<cellml_apiIReaction> rxn;
  rv = rxni->NextReaction(getter_AddRefs(rxn));
  while (rxn != NULL) {
    //parse the reaction;
    //reactions don't have names, so make up a new one for Antimony
    //LS DEBUG:  stopped coding here for now...
    //Variable* rxnvar = AddOrFindVariable(&cellmlname);
    rv = rxni->NextReaction(getter_AddRefs(rxn));
  }
  
  //Math
  nsCOMPtr<cellml_apiIMathContainer> mc(do_QueryInterface(component));
  nsCOMPtr<cellml_apiIMathList> mathlist;
  rv = mc->GetMath(getter_AddRefs(mathlist));
  nsCOMPtr<cellml_apiIMathMLElementIterator> mli;
  rv = mathlist->Iterate(getter_AddRefs(mli));
  nsCOMPtr<mathml_domIMathMLElement> mathel;
  rv = mli->Next(getter_AddRefs(mathel));
  while(mathel != NULL) {
    nsCOMPtr<domINodeList> nodes;
    rv = mathel->GetChildNodes(getter_AddRefs(nodes));
    PRUint32 length;
    rv = nodes->GetLength(&length);
    for (PRUint32 i=0; i<length; i++) {
      nsCOMPtr<domINode> node;
      rv = nodes->Item(i, getter_AddRefs(node));
      nsCOMPtr<mathml_domIMathMLApplyElement> input(do_QueryInterface(node, &rv));
      if (input != NULL) {
        nsCString indent;
        indent = "";
        nsCString cinfix;
        MathMLInputServices mmlis;
        rv = mmlis.MathMLToInputFormat(input, NULL, indent, cinfix);
        string infix;
        infix = cinfix.get();
        cout << infix << endl;
        string variable, equation;
        variable.assign(infix, 0, infix.find('='));
        equation.assign(infix, infix.find('=')+1, infix.size());
        variable = Trim(variable);
        equation = Trim(equation);
        string origeq = equation;

        //Remove all the '{id: ...}' bits
        size_t idpos;
        while ((idpos = equation.find("{id:")) != string::npos) {
          equation.erase(idpos, equation.find("}", idpos)-idpos+1);
        }
        //Remove '{unit: ...}' bits (for now)
        while ((idpos = equation.find("{unit:")) != string::npos) {
          equation.erase(idpos, equation.find("}", idpos)-idpos+1);
        }
        //Remove '{base: ...}' bits
        while ((idpos = equation.find("{base:")) != string::npos) {
          equation.erase(idpos, equation.find("}", idpos)-idpos+1);
        }
        //Remove '$'--it's apparantly the way that cellml notes variables that are also function names.
        // (We'll fix the name later ourselves with 'FixNames')
        while ((idpos = equation.find("$")) != string::npos) {
          equation.erase(idpos, 1);
        }
        while ((idpos = variable.find("$")) != string::npos) {
          variable.erase(idpos, 1);
        }

        //Handle piecewise equations
        if ((idpos = equation.find("piecewise")) != string::npos) {
          equation = CellMLPiecewiseToSBML(equation);
        }
        //Claim we can't handle definite integrals:
        if ((idpos = equation.find("definite_integral")) != string::npos) {
          g_registry.AddWarning("Unable to use the formula \"" + equation + "\" to define '" + variable + "' because Antimony does not handle definite integrals.");
          continue;
        }
        Formula* formula = g_registry.NewBlankFormula();
        setFormulaWithString(Trim(equation), formula, this);

        //Find out what variable we're assigning to, and how we're assigning to it.
        //Remove '{unit: ...}' bits (for now)
        while ((idpos = variable.find("{unit:")) != string::npos) {
          variable.erase(idpos, variable.find("}", idpos)-idpos+1);
        }
        vector<string> fullname;
        fullname.push_back(variable);
        Variable* var = GetVariable(fullname);
        if (var != NULL) {
          //Math is simple assignent rule
          if (var->SetAssignmentRule(formula)) {
            //Something went wrong
            //cout << "Unable to use the formula \"" << formula->ToDelimitedStringWithEllipses('.') << "\" (originally \"" << origeq << "\") to set the assignment rule for " << var->GetNameDelimitedBy('.') << ":  " << getLastError() << endl;
            string warning = "Unable to use the formula \"" + formula->ToDelimitedStringWithEllipses('.') + "\" to set the assignment rule for " + var->GetNameDelimitedBy('.') + ":  " + getLastError();
            g_registry.AddWarning(warning);
          }
        }
        else if (variable.find("d(") == 0 && variable.find(")/d(") != string::npos) {
          //Math is of the form dx/dy--if y is time, we can make this a rate rule.
          size_t timepos = variable.find(")/d(")+4;
          string maybetime;
          maybetime.assign(variable, timepos, variable.find(')', timepos)-timepos);
          cellmltext = ToNSString(maybetime);
          rv = varset->GetVariable(cellmltext, getter_AddRefs(cmlvar));
          if (cmlvar && !HasTimeUnits(cmlvar)) {
            rv = cmlvar->GetUnitsName(cellmltext);
            string warning = "The units of \"" + maybetime + "\" ('" + ToThinString(cellmltext.get()) + "') do not have 'seconds' as their base unit, so assuming this CellML model is trying to take the derivative of something with respect to some not-time element, we are not translating this derivative.";
            g_registry.AddWarning(warning);
            continue;
          }
          variable.assign(variable, 2, timepos-6);
          variable = Trim(variable);
          var = AddOrFindVariable(&variable);
          if (var->SetRateRule(formula)) {
            string warning = "Unable to use the formula \"" + formula->ToDelimitedStringWithEllipses('.') + "\" to set the rate rule for " + var->GetNameDelimitedBy('.') + ":  " + getLastError();
            g_registry.AddWarning(warning);
          }
        }
        else if (variable.find("del(") != string::npos) {
          //It's a partial differential equation.
          string warning = "Unable to translate an assignment to \"" + variable + "\" in the Antimony format because Antimony does not handle partial differential equations (though neither do existing CellML tools).";
          g_registry.AddWarning(warning);
        }
        else if (variable.find("selector(") != string::npos) {
          //It's vector or matrix math of some sort.
          string warning = "Unable to translate an assignment to \"" + variable + "\" in the Antimony format because Antimony does not handle vector or matrix algebra.";
          g_registry.AddWarning(warning);
        }
        else if (IsReal(variable)) {
          //It's some sort of algebraic rule.
          string warning = "Unable to translate the equation \"" + variable + " = " + equation + "\" because Antimony does not handle algebraic rules.";
          g_registry.AddWarning(warning);
        }
        else {
          //Unable to determine what kind of math we're talking about.
          string warning = "Unable to figure out how to translate an assignment to \"" + variable + "\" in the Antimony format.  This variable may have been left undefined in that component.";
          g_registry.AddWarning(warning);
        }

      }
      else {
        //string warning = "Child node " +  i + " of 'math' element not an 'apply' element";
      }
    }
    rv = mli->Next(getter_AddRefs(mathel));
  }

  //Containers (?)

  //And finally, fix names.
  FixNames();
}

void Module::SetCellMLChildrenAsSubmodules(nsCOMPtr<cellml_apiICellMLComponent> component) {
  //Iterate over 'encapsulation' children and make them submodules.
  nsString cellmltext;
  string cellmlname;
  nsresult rv;
  nsCOMPtr<cellml_apiICellMLComponentSet> children;
  rv = component->GetEncapsulationChildren(getter_AddRefs(children));
  nsCOMPtr<cellml_apiICellMLComponentIterator> childi;
  nsCOMPtr<cellml_apiICellMLComponent> child;
  
  rv = children->IterateComponents(getter_AddRefs(childi));
  rv = childi->NextComponent(getter_AddRefs(child));
  while (child != NULL) {
    cellmlname = GetNameAccordingToEncapsulationParent(child, m_cellmlmodel);
    //rv = child->GetName(cellmltext);
    //cellmlname = ToThinString(cellmltext.get());
    //FixName(cellmlname);
    string cellmlmodname = GetModuleNameFrom(child);
    Module* submod = g_registry.GetModule(cellmlmodname);
    submod->SetCellMLChildrenAsSubmodules(child); //Recursive!  This is so the submodels are all set up before we copy them.
    vector<string> fullname;
    fullname.push_back(cellmlname);
    Variable* foundvar = GetVariable(fullname);
    if (foundvar != NULL && !(foundvar->GetType()==varModule && m_childrenadded)) {
      cellmlname = cellmlname + "_mod";
    }
    //Save the name, since it's not obvious whether the "_mod" was added or not.
    nsCString compmodid;
    nsCOMPtr<IWrappedPCM> cw(do_QueryInterface(child));
    rv = cw->GetObjid(compmodid);
    g_registry.m_cellmlnames.insert(make_pair(compmodid.get(), cellmlname)); //Even if we've already added this submodule, each time it's imported, the submodule gets its own component ID, and they all need to go in here.
    if (!m_childrenadded) {
      Variable* var = AddOrFindVariable(&cellmlname);
      if(var->SetModule(&cellmlmodname)) {
        assert(false);
        return;
      }
    }
    rv = childi->NextComponent(getter_AddRefs(child));
  }
  m_childrenadded = true; //Since this is recursive, we may call some multiply-imported submodules multiple times otherwise.
}


const nsCOMPtr<cellml_apiIModel> Module::GetCellMLModel()
{
  if (m_cellmlmodel==NULL) {
    CreateCellMLModel();
  }
  else {
    nsString cellmltext;
    m_cellmlmodel->GetName(cellmltext);
    if (ToThinString(cellmltext.get()) != m_modulename) {
      CreateCellMLModel();
    }
  }
  return m_cellmlmodel;
}

void Module::ReloadSubmodelVariables(const string& modulename)
{
  for (size_t var=0; var<m_variables.size(); var++) {
    Variable* variable = m_variables[var];
    if (variable->GetType() == varModule) {
      Module* modcopy = variable->GetModule();
      Module* submod = g_registry.GetModule(modcopy->GetModuleName());
      submod->ReloadSubmodelVariables(submod->GetModuleName());
      modcopy->ReloadSubmodelVariables(modulename);
      modcopy->ResyncVariablesWith(submod, modulename, variable->GetName());
    }
  }
}

void Module::ResyncVariablesWith(const Module* twin, string modulename, vector<string> varname)
{
  assert (m_variables.size() <= twin->m_variables.size());
  for (size_t var = m_variables.size(); var<twin->m_variables.size(); var++) {
    //New variables
    Variable* newsubvar = twin->m_variables[var];
    Variable* newvar = new Variable(*newsubvar);
    for (size_t name=varname.size(); name>0; name--) {
      newvar->SetNewTopName(modulename, varname[name-1]);
    }
    //cout << "new subvar: " << ToStringFromVecDelimitedBy(newvar->GetName(), '.') << " for module " << ToStringFromVecDelimitedBy(m_variablename, '.') << endl;
    assert(newvar->GetType() != varModule);
    m_variables.push_back(newvar);
    StoreVariable(newvar);
  }
}

void Module::ReloadSubmodelConnections(Module* syncmod)
{
  for (size_t var=0; var<m_variables.size(); var++) {
    Variable* variable = m_variables[var];
    if (variable->GetType() == varModule) {
      Module* modcopy = variable->GetModule();
      Module* submod = g_registry.GetModule(modcopy->GetModuleName());
      assert(modcopy != submod);
      submod->ReloadSubmodelConnections(submod);
      modcopy->ReloadSubmodelConnections(syncmod);
      //cout << "reloading connections from " << submod->m_modulename << " to sync with " << variable->GetNameDelimitedBy('.') << endl;
      m_varmap.insert(submod->m_varmap.begin(), submod->m_varmap.end());
      for (size_t sync = modcopy->m_synchronized.size(); sync<submod->m_synchronized.size(); sync++) {
        //There are synchronizations in the ur-module that we didn't get when we copied it.
        vector<string> var1name = submod->m_synchronized[sync].first;
        vector<string> var2name = submod->m_synchronized[sync].second;
        for (size_t vn=variable->GetName().size(); vn>0; vn--) {
          var1name.insert(var1name.begin(), variable->GetName()[vn-1]);
          var2name.insert(var2name.begin(), variable->GetName()[vn-1]);
        }
        Variable* var1 = syncmod->GetVariable(var1name);
        Variable* var2 = syncmod->GetVariable(var2name);
        size_t numsynced = syncmod->m_synchronized.size();
        assert(var1 != NULL && var2 != NULL);
        if (var1->GetIsEquivalentTo(var2)) {
          //Already synchronized.
          continue;
        }
        if (!var1->Synchronize(var2)) {
          //This adds the synchronization to the local list instead of the submodel's list.  So, move it!
          assert(g_registry.GetModule(var1->GetNamespace()) == syncmod);
          assert(numsynced == syncmod->m_synchronized.size()-1);
          assert(syncmod->m_synchronized.size() > 0);
          pair<vector<string>, vector<string> > newsync = syncmod->m_synchronized[syncmod->m_synchronized.size()-1];
          syncmod->m_synchronized.pop_back();
          modcopy->m_synchronized.push_back(newsync);
          /*
          cout << ToStringFromVecDelimitedBy(var1name, '.') << " to "
               << ToStringFromVecDelimitedBy(var2name, '.') << ": " << endl;
          cout << ToStringFromVecDelimitedBy(var1->GetName(), '.') << " and "
               << ToStringFromVecDelimitedBy(var2->GetName(), '.') << " now synchronized";
          if (!var1->GetFormula()->IsEmpty()) {
            cout << ": " << var1->GetFormula()->ToDelimitedStringWithEllipses('.') << endl;
          }
          else {
            cout << "." << endl;
          }
          */
        }
        else {
          g_registry.AddWarning("In module '" + m_modulename + "', the variables " + var1->GetNameDelimitedBy('.') + " and " + var2->GetNameDelimitedBy('.') + " were unable to be set as equivalent:  " + g_registry.GetError());
          //This is tricky, because now we have to remove the main module's synchronization that set this up in the first place.
          vector<string> sync1name, sync2name;
          bool foundorig = false;
          if (var1->IsPointer()) {
            cout << "var1 is pointer: " << var1->GetNameDelimitedBy('.') << endl;
            sync1name = var1->GetName();
            sync2name = var1->GetPointerName();
            cout << "Looking for sync for " << ToStringFromVecDelimitedBy(sync1name, '.') << " and " << ToStringFromVecDelimitedBy(sync2name, '.') << endl;
            for (vector<pair<vector<string>, vector<string> > >::iterator sync = m_synchronized.begin();
                 sync != m_synchronized.end(); sync++) {
              cout << "Found " << ToStringFromVecDelimitedBy(sync->first, '.') << " and " << ToStringFromVecDelimitedBy(sync->second, '.') << endl;
              if ((sync1name == sync->first && sync2name == sync->second) ||
                  (sync2name == sync->first && sync1name == sync->second)) {
                //cout << "Removing sync from " << m_modulename << ": " << ToStringFromVecDelimitedBy(sync1name, '.') << " with " << ToStringFromVecDelimitedBy(sync2name, '.') << endl;
                m_synchronized.erase(sync);
                foundorig = true;
                //cout << "Found original!" << endl;
                break;
              }
            }
          }
          if (!foundorig && var2->IsPointer()) {
            cout << "var2 is pointer: " << var2->GetNameDelimitedBy('.') << endl;
            sync1name = var2->GetName();
            sync2name = var2->GetPointerName();
            for (vector<pair<vector<string>, vector<string> > >::iterator sync = m_synchronized.begin();
                 sync != m_synchronized.end(); sync++) {
              if ((sync1name == sync->first && sync2name == sync->second) ||
                  (sync2name == sync->first && sync1name == sync->second)) {
                //cout << "Removing sync from " << m_modulename << ": " << ToStringFromVecDelimitedBy(sync1name, '.') << " with " << ToStringFromVecDelimitedBy(sync2name, '.') << endl;
                m_synchronized.erase(sync);
                foundorig = true;
                //cout << "Found original!" << endl;
                break;
              }
            }
          }
          
          assert(foundorig);
        }
      }
    }
  }
  //Finalize(); //Re-finalize to update m_origvars.
}

//
//----------Creating CellML Models---------------
//

void Module::CreateCellMLModel()
{
  nsresult rv;
  nsString cellmltext;
  wstring cellmlwstring;

  if (m_cellmlmodel != NULL) {
    if (m_cellmlcomponent != NULL) {
      m_cellmlcomponent = NULL;
    }
  }

  nsCOMPtr<cellml_apiICellMLBootstrap> boot(do_GetService(CELLML_BOOTSTRAP_CONTRACTID, &rv));
  NS_ENSURE_SUCCESS(rv, );
  rv = boot->CreateModel(NS_LITERAL_STRING("1.1"), getter_AddRefs(m_cellmlmodel));
  NS_ENSURE_SUCCESS(rv, );
  nsCOMPtr<domIDOMImplementation> domi;
  rv = boot->GetDomImplementation(getter_AddRefs(domi));
  assert(domi != NULL);
  nsCOMPtr<domIDocument> doc;
  nsCOMPtr<domIDocumentType> dt;
  nsString model = ToNSString("model");
  nsString none =  ToNSString("");
  rv = domi->CreateDocumentType(model, none, none, getter_AddRefs(dt));
  rv = domi->CreateDocument(none, model, dt, getter_AddRefs(doc));

  cellmltext = ToNSString(m_modulename);
  m_cellmlmodel->SetName(cellmltext);
  //Create units

  //Create a component for all local variables
  AddCellMLComponentsTo(m_cellmlmodel, this);

  //Create encapsulation relationships
  AddEncapsulationTo(m_cellmlmodel);

  //Create all connections
  AddConnectionsTo(m_cellmlmodel, this);

  //Add in the Math
  for (map<Variable*, vector<Variable*> >::iterator mapiter = m_syncedvars.begin();
       mapiter != m_syncedvars.end(); mapiter++) {
    AssignMathOnceFor(mapiter->second, doc);
  }

  //Add in the ODEs
  AddODEsTo(m_cellmlmodel, this);


}

void Module::AddCellMLComponentsTo(nsCOMPtr<cellml_apiIModel> model, Module* topmod)
{
  model->AddElement(GetCellMLComponent(topmod));
  for (size_t var=0; var<m_variables.size(); var++) { 
    Variable* variable = m_variables[var];
    if (variable->GetType() == varModule) {
      variable->GetModule()->AddCellMLComponentsTo(model, topmod);
    }
  }
}

nsCOMPtr<cellml_apiICellMLComponent> Module::GetCellMLComponent(Module* topmod)
{
  if (m_cellmlcomponent == NULL) {
    CreateCellMLComponent(topmod);
  }
  return m_cellmlcomponent;
}

void Module::CreateCellMLComponent(Module* topmod)
{
  nsresult rv;
  //Establish a unique name in CellML, which is a different namespace than Antimony (which can distinguish X.y from Z.y--in CelML, we can no longer call both 'y', and it's awkward to call everything 'X_y', etc.)
  int nindex = static_cast<int>(m_variablename.size())-1;
  string name = m_modulename; //The top component needs a name here, too.
  if (nindex > 0) {
    name = m_variablename[nindex];
    while (topmod->InUnique(name)) {
      assert(nindex>0);
      if (nindex==0) {
        //error  LS DEBUG
      }
      nindex--;
      name = m_variablename[nindex] + "_" + name;
    }
  }
  topmod->AddUnique(m_variablename, name);
  rv = topmod->m_cellmlmodel->CreateComponent(getter_AddRefs(m_cellmlcomponent));
  rv = m_cellmlcomponent->SetName(ToNSString(name));
  m_cellmlmodel = topmod->m_cellmlmodel;

  map<Variable*, vector<Variable*> >::iterator mapiter;
  
  //Create into m_cellmlcomponent.
  for (size_t var=0; var<m_variables.size(); var++) { 
    Variable* variable = m_variables[var];
    bool used = false;
    switch(variable->GetType()) {
    case varSpeciesUndef:
    case varFormulaUndef:
    case varDNA:
    case varFormulaOperator:
    case varReactionGene:
    case varInteraction:
    case varReactionUndef:
    case varCompartment:
    case varUndefined:
      //The above all become CellML variables:
      AddVariableToCellML(variable, topmod->m_cellmlmodel);
      used = true;
      break;
    case varEvent:
    case varModule:
    case varStrand:
      //These all have no CellML equivalent (except Module, which is taken care of separately).
      break;
    }
    if (used) {
      //Create a list of synchronized variables
      Variable* samevar = variable->GetSameVariable();
      mapiter = topmod->m_syncedvars.find(samevar);
      if (mapiter == topmod->m_syncedvars.end()) {
        vector<Variable*> varlist;
        varlist.push_back(variable);
        topmod->m_syncedvars.insert(make_pair(samevar, varlist));
      }
      else {
        mapiter->second.push_back(variable);
      }
    }
  }
}

void Module::AddVariableToCellML(Variable* variable, nsCOMPtr<cellml_apiIModel> model)
{
  vector<string> varname = variable->GetName();
  assert(varname.size()>0);
  nsCOMPtr<cellml_apiICellMLVariable> cmlvar = AddVariableToCellML(varname[varname.size()-1], model);
  variable->SetCellMLVariable(cmlvar);
}

nsCOMPtr<cellml_apiICellMLVariable> Module::AddVariableToCellML(string varname, nsCOMPtr<cellml_apiIModel> model)
{
  assert(m_cellmlcomponent != NULL);
  nsresult rv;
  nsCOMPtr<cellml_apiICellMLVariableSet> cmlvarset;
  rv = m_cellmlcomponent->GetVariables(getter_AddRefs(cmlvarset));
  nsCOMPtr<cellml_apiICellMLVariable> cmlvar;
  nsString cmlvarst = ToNSString(varname);
  rv = cmlvarset->GetVariable(cmlvarst, getter_AddRefs(cmlvar));
  if (cmlvar != NULL) {
    //Already exists!
    return cmlvar;
  }
  rv = model->CreateCellMLVariable(getter_AddRefs(cmlvar));
  rv = m_cellmlcomponent->AddElement(cmlvar);
  rv = cmlvar->SetName(ToNSString(varname));
  rv = cmlvar->SetPublicInterface(1);
  rv = cmlvar->SetPrivateInterface(0); //LS DEBUG: revisit this later
  return cmlvar;
}

void Module::AssignMathOnceFor(vector<Variable*> varlist, nsCOMPtr<domIDocument> doc)
{
  nsresult rv;
  Variable* finalvar = varlist[0];
  const Formula* ia = finalvar->GetInitialAssignment();
  const Formula* ar = finalvar->GetAssignmentRuleOrKineticLaw();
  const Formula* rr = finalvar->GetRateRule();
  //All three may have been initially defined in different modules
  if (!ia->IsEmpty()) {
    if (ia->IsDouble()) {
      Variable* targetvar = WhichFirstDefined(varlist, formulaINITIAL);
      nsCOMPtr<cellml_apiICellMLVariable> cmlvar = targetvar->GetCellMLVariable();
      rv = cmlvar->SetInitialValue(ToNSString(ia->ToDelimitedStringWithEllipses('_')));
      cout << "Successfully set initial value for " << targetvar->GetNameDelimitedBy('.') << endl;
    }
    else {
      //Can't do this in CellML.
      cout << "Can't set initial assignment for " << varlist[0]->GetNameDelimitedBy('.') << endl;
    }
  }
  if (!ar->IsEmpty()) {
    Variable* targetvar = WhichFirstDefined(varlist, formulaASSIGNMENT);
    const Variable* origtarget = targetvar->GetOriginal();
    vector<string> varname = origtarget->GetName();
    assert(varname.size()==1);
    string formula = origtarget->GetAssignmentRuleOrKineticLaw()->ToDelimitedStringWithStrands('_', origtarget->GetStrandVars());
    formula = varname[varname.size()-1] + " = " + formula;
    if (AddCellMLMathTo(formula, targetvar, doc)) {
      cout << "Successfully found assignment rule for " << targetvar->GetNameDelimitedBy('.') << " (" << formula << ")" << endl;
    }
    else {
      cout << "Unable to translate \"" << formula << "\" to CellML's MathML for the assignment rule." << endl;
    }
  }
  if (!rr->IsEmpty()) {
    Variable* targetvar = WhichFirstDefined(varlist, formulaRATE);
    AddTimeFor(targetvar->GetCellMLVariable());
    const Variable* origtarget = targetvar->GetOriginal();
    vector<string> varname = origtarget->GetName();
    assert(varname.size()==1);
    string formula = origtarget->GetRateRule()->ToDelimitedStringWithStrands('_', origtarget->GetStrandVars());
    formula = "d(" + varname[varname.size()-1] + ")/d(time) = " + formula;
    if (AddCellMLMathTo(formula, targetvar, doc)) {
      cout << "Successfully found rate rule for " << targetvar->GetNameDelimitedBy('.') << " (" << formula << ")" << endl;
    }
    else {
      cout << "Unable to translate \"" << formula << "\" to CellML's MathML for the rate rule." << endl;
    }
  }
}

bool Module::AddCellMLMathTo(string formula, Variable* targetvar, nsCOMPtr<domIDocument> doc)
{
  nsresult rv;
  nsCOMPtr<cellml_apiICellMLVariable> cmlvar = targetvar->GetCellMLVariable();
  nsCOMPtr<cellml_apiICellMLComponent> cmlcomp = GetCellMLComponentOf(cmlvar);
  nsCOMPtr<cellml_apiIMathContainer> upcast_cmlcomp;
  rv = cmlcomp->QueryInterface(upcast_cmlcomp->GetIID(), getter_AddRefs(upcast_cmlcomp));
  nsCOMPtr<domIElement> mathml;
  nsCString formstring = ToNSCString(formula);
  MathMLInputServices mmlis;
  rv = mmlis.InputFormatToMathML(doc, formstring, getter_AddRefs(mathml));
  NS_ENSURE_SUCCESS(rv, false);
  nsCOMPtr<mathml_domIMathMLElement> recast_mathml;
  rv = mathml->QueryInterface(recast_mathml->GetIID(), getter_AddRefs(recast_mathml));
  NS_ENSURE_SUCCESS(rv, false);
  rv = upcast_cmlcomp->AddMath(recast_mathml);
  NS_ENSURE_SUCCESS(rv, false);
  return true;
}

void Module::AddTimeFor(nsCOMPtr<cellml_apiICellMLVariable> cmlvar)
{
  nsCOMPtr<cellml_apiICellMLComponent> cmlcomp = GetCellMLComponentOf(cmlvar);
  AddTimeTo(cmlcomp);
}

nsCOMPtr<cellml_apiICellMLVariable> Module::AddTimeTo(nsCOMPtr<cellml_apiICellMLComponent> cmlcomp)
{
  nsresult rv;
  nsCOMPtr<cellml_apiICellMLVariableSet> cmlvarset;
  rv = cmlcomp->GetVariables(getter_AddRefs(cmlvarset));
  nsCOMPtr<cellml_apiICellMLVariable> time;
  nsString timest = ToNSString("time");
  rv = cmlvarset->GetVariable(timest, getter_AddRefs(time));
  if (time != NULL) {
    //Already exists!
    return time;
  }
  assert(m_cellmlmodel != NULL);
  rv = m_cellmlmodel->CreateCellMLVariable(getter_AddRefs(time));
  rv = cmlcomp->AddElement(time);
  rv = time->SetName(timest);
  rv = time->SetPublicInterface(1);
  rv = time->SetPrivateInterface(0); //LS DEBUG: revisit this later

  nsCOMPtr<cellml_apiICellMLComponent> parent;
  rv = cmlcomp->GetEncapsulationParent(getter_AddRefs(parent));
  if (parent != NULL) {
    nsCOMPtr<cellml_apiICellMLVariable> ptime = AddTimeTo(parent);
    AddOneConnection(m_cellmlmodel, time, ptime);
  }
  return time;
}

Variable* Module::WhichFirstDefined(vector<Variable*> varlist, formula_type ftype)
{
  assert(varlist.size()>0);
  Variable* candidate = NULL;
  for (size_t var=0; var<varlist.size(); var++) {
    if (varlist[var]->StillMatchesOriginal(ftype)) {
      if (candidate == NULL || candidate->GetName().size() < varlist[var]->GetName().size()) {
        candidate = varlist[var];
      }
    }
  }
  if (candidate == NULL) {
    //The top module changed a variable's definition without giving it a local name.
    //Um, help? LS DEBUG
    //At least point to the version that everyone is synced to:
    candidate = varlist[0]->GetSameVariable();
    cout << "Help!  Nothing found for the variables synced to " << candidate->GetNameDelimitedBy('.');
  }
  return candidate;
}

bool Module::InUnique(string name)
{
  return (m_uniquenames.find(name) != m_uniquenames.end());
}

void Module::AddUnique(vector<string> fullname, string name)
{
  m_uniquenames.insert(name);
  m_cellmlnames.insert(make_pair(fullname, name));
}

string Module::GetCellMLNameOf(vector<string> fullname)
{
  assert(m_cellmlnames.find(fullname) != m_cellmlnames.end());
  return m_cellmlnames.find(fullname)->second;
}

void Module::AddEncapsulationTo(nsCOMPtr<cellml_apiIModel> model)
{
  nsresult rv;
  nsCOMPtr<cellml_apiIGroup> group;
  rv = model->CreateGroup(getter_AddRefs(group));
  rv = model->AddElement(group);
  nsCOMPtr<cellml_apiIRelationshipRef> relref;
  rv = model->CreateRelationshipRef(getter_AddRefs(relref));
  rv = group->AddElement(relref);
  rv = relref->SetRelationshipName(ToNSString(""), ToNSString("encapsulation"));
  vector<string> blank;
  nsCOMPtr<cellml_apiIComponentRef> cr = GetComponentRef(m_cellmlmodel, GetCellMLNameOf(blank));
  rv = group->AddElement(cr);
}

nsCOMPtr<cellml_apiIComponentRef> Module::GetComponentRef(nsCOMPtr<cellml_apiIModel> model, string cmlname)
{
  nsresult rv;
  nsCOMPtr<cellml_apiIComponentRef> cr;
  rv = model->CreateComponentRef(getter_AddRefs(cr));
  rv = cr->SetComponentName(ToNSString(cmlname));
  for (size_t var=0; var<m_variables.size(); var++) {
    if (m_variables[var]->GetType() == varModule) {
      string subvarcmlname = GetCellMLNameOf(m_variables[var]->GetName());
      nsCOMPtr<cellml_apiIComponentRef> subcr = m_variables[var]->GetModule()->GetComponentRef(model, subvarcmlname);
      rv = cr->AddElement(subcr);
    }
  }
  return cr;
}

void Module::AddConnectionsTo(nsCOMPtr<cellml_apiIModel> model, Module* topmod)
{
  for (size_t var=0; var<m_variables.size(); var++) {
    Variable* var1 = m_variables[var];
    if (var1->IsPointer()) {
      AddOneConnection(model, var1, topmod);
    }
    else if (var1->GetType()==varModule) {
      var1->GetModule()->AddConnectionsTo(model, topmod);
    }
  }
}

void Module::AddOneConnection(nsCOMPtr<cellml_apiIModel> model, Variable* var, Module* topmod)
{
  nsresult rv;
  assert(var->IsPointer());
  nsCOMPtr<cellml_apiIConnection> connection;
  rv = model->CreateConnection(getter_AddRefs(connection));

  //Variables
  vector<string> var1name = var->GetName();
  vector<string> var2name = var->GetSameVariable()->GetName();
  string cellmlname1 = var1name[var1name.size()-1];
  string cellmlname2 = var2name[var2name.size()-1];
  nsCOMPtr<cellml_apiIMapVariables> mapvars;
  rv = model->CreateMapVariables(getter_AddRefs(mapvars));
  rv = connection->AddElement(mapvars);
  rv = mapvars->SetFirstVariableName(ToNSString(cellmlname1));
  rv = mapvars->SetSecondVariableName(ToNSString(cellmlname2));

  //Components
  vector<string> mod1name = var1name;
  mod1name.pop_back();
  vector<string> mod2name = var2name;
  mod2name.pop_back();
  string cellmlcomp1 = topmod->GetCellMLNameOf(mod1name);
  string cellmlcomp2 = topmod->GetCellMLNameOf(mod2name);
  nsCOMPtr<cellml_apiIMapComponents> compmap;
  rv = connection->GetComponentMapping(getter_AddRefs(compmap));
  rv = compmap->SetFirstComponentName(ToNSString(cellmlcomp1));
  rv = compmap->SetSecondComponentName(ToNSString(cellmlcomp2));
  rv = model->AddElement(connection);
}

void Module::AddOneConnection(nsCOMPtr<cellml_apiIModel> model, nsCOMPtr<cellml_apiICellMLVariable> var1, nsCOMPtr<cellml_apiICellMLVariable> var2)
{
  nsresult rv;
  nsCOMPtr<cellml_apiIConnection> connection;
  rv = model->CreateConnection(getter_AddRefs(connection));
  rv = model->AddElement(connection);
  nsCOMPtr<cellml_apiIMapVariables> mapvars;
  rv = model->CreateMapVariables(getter_AddRefs(mapvars));
  rv = connection->AddElement(mapvars);
  rv = mapvars->SetFirstVariable(var1);
  rv = mapvars->SetSecondVariable(var2);
}

void Module::AddODEsTo(nsCOMPtr<cellml_apiIModel> model, Module* topmod)
{
  set<Variable*> species;
  set<Variable*> reactions;
  GetAllSpeciesAndReactions(species, reactions);
  set<Variable*>::iterator speciter, rxniter;
  for (speciter=species.begin(); speciter != species.end(); speciter++) {
    vector<string> commonmod;
    Formula form;
    set<Variable*> involvedrxns;
    for (rxniter=reactions.begin(); rxniter != reactions.end(); rxniter++) {
      double stoich = (*rxniter)->GetReaction()->GetStoichiometryFor(*speciter);
      if (stoich != 0) {
        if (stoich < 0) {
          form.AddMathThing('-');
        }
        else if (!form.IsEmpty()) {
          form.AddMathThing('+');
        }
        if (stoich != 1.0 && stoich != -1.0) {
          form.AddNum(stoich);
          form.AddMathThing('*');
        }
        form.AddVariable(*rxniter);
        involvedrxns.insert(*rxniter);
      }
    }
    if (involvedrxns.size()==0) continue; //The species was in no reactions.
    set<Variable*> contains;
    Module* ratemod = topmod->BestModuleToAdd(involvedrxns, contains);
    ratemod->AddRateRuleInvolving(*speciter, form, involvedrxns);
  }
}

void Module::GetAllSpeciesAndReactions(set<Variable*>& species, set<Variable*>& reactions)
{
  for (size_t var=0; var<m_variables.size(); var++) {
    if (IsSpecies(m_variables[var]->GetType()) && !(m_variables[var]->GetIsConst())) {
      species.insert(m_variables[var]->GetSameVariable());
    }
    else if (IsReaction(m_variables[var]->GetType())) {
      reactions.insert(m_variables[var]->GetSameVariable());
    }
    else if (m_variables[var]->GetType()==varModule) {
      m_variables[var]->GetModule()->GetAllSpeciesAndReactions(species, reactions);
    }
  }
}

Module* Module::BestModuleToAdd(set<Variable*> involvedrxns, set<Variable*>& contains )
{
  for (size_t var=0; var<m_variables.size(); var++) {
    var_type vtype = m_variables[var]->GetType();
    if (vtype==varModule) {
      set<Variable*> partcontains;
      Module* testmod = m_variables[var]->GetModule()->BestModuleToAdd(involvedrxns, partcontains);
      assert(testmod == NULL || (partcontains.size() == involvedrxns.size()));
        set<Variable*>::iterator contit;
      for (contit = partcontains.begin(); contit != partcontains.end(); contit++) {
        contains.insert(*contit);
      }
      if (testmod != NULL) {
        return testmod; //Might be some parallel sub-module that *also* contains all rxns, but who cares.
      }
    }
    else if (IsReaction(vtype)) {
      Variable* sourcerxn = m_variables[var]->GetSameVariable();
      if (involvedrxns.find(sourcerxn) != involvedrxns.end()) {
        contains.insert(sourcerxn);
      }
    }
  }
  if (involvedrxns.size() == contains.size()) return this;
  return NULL;
}

void Module::AddRateRuleInvolving(Variable* species, Formula form, set<Variable*> involvedrxns)
{
  nsCOMPtr<cellml_apiICellMLVariable> subvar;
  string localname = FindOrCreateLocalVersionOf(species, subvar);
  for (set<Variable*>::iterator involvedit=involvedrxns.begin(); involvedit != involvedrxns.end(); involvedit++) {
    string localrxn = FindOrCreateLocalVersionOf(*involvedit, subvar);
    form.UseInstead(localrxn, *involvedit);
  }
  string infix = "d(" + localname + ")/d(time) = " + form.ToCellML();
  AddTimeTo(m_cellmlcomponent);
  cout << "If we could convert infix to CellML yet, we would get: " << infix << endl;
}

string Module::FindOrCreateLocalVersionOf(Variable* variable, nsCOMPtr<cellml_apiICellMLVariable>& localvar)
{
  for (size_t var=0; var<m_variables.size(); var++) {
    if (variable->GetSameVariable() == m_variables[var]->GetSameVariable()) {
      vector<string> varname = m_variables[var]->GetName();
      //assert(varname.size()==1); //Not true!  Names have *this* module's name in the front.
      localvar = m_variables[var]->GetCellMLVariable();
      return (varname[varname.size()-1]);
    }
  }
  for (size_t var=0; var<m_variables.size(); var++) {
    if (m_variables[var]->GetType()==varModule) {
      nsCOMPtr<cellml_apiICellMLVariable> subvar;
      string foundvar = m_variables[var]->GetModule()->FindOrCreateLocalVersionOf(variable, subvar);
      if (foundvar != "") {
        //The variable was indeed in this list.  Create a local copy (in CellML) and sync it.
        vector<string> varname;
        varname.push_back(foundvar);
        while (GetVariable(varname) != NULL) {
          foundvar = m_variables[var]->GetModule()->GetModuleName() + "_" + foundvar;
          varname[0] = foundvar;
        }
        localvar = AddVariableToCellML(foundvar, m_cellmlmodel);
        AddOneConnection(m_cellmlmodel, localvar, subvar);
        return foundvar;
      }
    }
  }
  return "";
}

#endif
