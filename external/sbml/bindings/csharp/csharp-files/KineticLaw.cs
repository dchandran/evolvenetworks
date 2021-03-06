/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

namespace libsbml {

using System;
using System.Runtime.InteropServices;

public class KineticLaw : SBase {
	private HandleRef swigCPtr;
	
	internal KineticLaw(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.KineticLawUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.KineticLawUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(KineticLaw obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (KineticLaw obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~KineticLaw() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_KineticLaw(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public KineticLaw(long level, long version) : this(libsbmlPINVOKE.new_KineticLaw__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public KineticLaw(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_KineticLaw__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public KineticLaw(KineticLaw orig) : this(libsbmlPINVOKE.new_KineticLaw__SWIG_2(KineticLaw.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public new KineticLaw clone() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_clone(swigCPtr);
    KineticLaw ret = (cPtr == IntPtr.Zero) ? null : new KineticLaw(cPtr, true);
    return ret;
  }

  public string getFormula() {
    string ret = libsbmlPINVOKE.KineticLaw_getFormula(swigCPtr);
    return ret;
  }

  public ASTNode getMath() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getMath(swigCPtr);
    ASTNode ret = (cPtr == IntPtr.Zero) ? null : new ASTNode(cPtr, false);
    return ret;
  }

  public string getTimeUnits() {
    string ret = libsbmlPINVOKE.KineticLaw_getTimeUnits(swigCPtr);
    return ret;
  }

  public string getSubstanceUnits() {
    string ret = libsbmlPINVOKE.KineticLaw_getSubstanceUnits(swigCPtr);
    return ret;
  }

  public bool isSetFormula() {
    bool ret = libsbmlPINVOKE.KineticLaw_isSetFormula(swigCPtr);
    return ret;
  }

  public bool isSetMath() {
    bool ret = libsbmlPINVOKE.KineticLaw_isSetMath(swigCPtr);
    return ret;
  }

  public bool isSetTimeUnits() {
    bool ret = libsbmlPINVOKE.KineticLaw_isSetTimeUnits(swigCPtr);
    return ret;
  }

  public bool isSetSubstanceUnits() {
    bool ret = libsbmlPINVOKE.KineticLaw_isSetSubstanceUnits(swigCPtr);
    return ret;
  }

  public int setFormula(string formula) {
    int ret = libsbmlPINVOKE.KineticLaw_setFormula(swigCPtr, formula);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setMath(ASTNode math) {
    int ret = libsbmlPINVOKE.KineticLaw_setMath(swigCPtr, ASTNode.getCPtr(math));
    return ret;
  }

  public int setTimeUnits(string sid) {
    int ret = libsbmlPINVOKE.KineticLaw_setTimeUnits(swigCPtr, sid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setSubstanceUnits(string sid) {
    int ret = libsbmlPINVOKE.KineticLaw_setSubstanceUnits(swigCPtr, sid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int unsetTimeUnits() {
    int ret = libsbmlPINVOKE.KineticLaw_unsetTimeUnits(swigCPtr);
    return ret;
  }

  public int unsetSubstanceUnits() {
    int ret = libsbmlPINVOKE.KineticLaw_unsetSubstanceUnits(swigCPtr);
    return ret;
  }

  public int addParameter(Parameter p) {
    int ret = libsbmlPINVOKE.KineticLaw_addParameter(swigCPtr, Parameter.getCPtr(p));
    return ret;
  }

  public int addLocalParameter(LocalParameter p) {
    int ret = libsbmlPINVOKE.KineticLaw_addLocalParameter(swigCPtr, LocalParameter.getCPtr(p));
    return ret;
  }

  public Parameter createParameter() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_createParameter(swigCPtr);
    Parameter ret = (cPtr == IntPtr.Zero) ? null : new Parameter(cPtr, false);
    return ret;
  }

  public LocalParameter createLocalParameter() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_createLocalParameter(swigCPtr);
    LocalParameter ret = (cPtr == IntPtr.Zero) ? null : new LocalParameter(cPtr, false);
    return ret;
  }

  public ListOfParameters getListOfParameters() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getListOfParameters__SWIG_0(swigCPtr);
    ListOfParameters ret = (cPtr == IntPtr.Zero) ? null : new ListOfParameters(cPtr, false);
    return ret;
  }

  public ListOfLocalParameters getListOfLocalParameters() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getListOfLocalParameters__SWIG_0(swigCPtr);
    ListOfLocalParameters ret = (cPtr == IntPtr.Zero) ? null : new ListOfLocalParameters(cPtr, false);
    return ret;
  }

  public Parameter getParameter(long n) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getParameter__SWIG_0(swigCPtr, n);
    Parameter ret = (cPtr == IntPtr.Zero) ? null : new Parameter(cPtr, false);
    return ret;
  }

  public LocalParameter getLocalParameter(long n) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getLocalParameter__SWIG_0(swigCPtr, n);
    LocalParameter ret = (cPtr == IntPtr.Zero) ? null : new LocalParameter(cPtr, false);
    return ret;
  }

  public Parameter getParameter(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getParameter__SWIG_2(swigCPtr, sid);
    Parameter ret = (cPtr == IntPtr.Zero) ? null : new Parameter(cPtr, false);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public LocalParameter getLocalParameter(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getLocalParameter__SWIG_2(swigCPtr, sid);
    LocalParameter ret = (cPtr == IntPtr.Zero) ? null : new LocalParameter(cPtr, false);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public long getNumParameters() { return (long)libsbmlPINVOKE.KineticLaw_getNumParameters(swigCPtr); }

  public long getNumLocalParameters() { return (long)libsbmlPINVOKE.KineticLaw_getNumLocalParameters(swigCPtr); }

  public UnitDefinition getDerivedUnitDefinition() {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_getDerivedUnitDefinition__SWIG_0(swigCPtr);
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, false);
    return ret;
  }

  public bool containsUndeclaredUnits() {
    bool ret = libsbmlPINVOKE.KineticLaw_containsUndeclaredUnits__SWIG_0(swigCPtr);
    return ret;
  }

  public Parameter removeParameter(long n) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_removeParameter__SWIG_0(swigCPtr, n);
    Parameter ret = (cPtr == IntPtr.Zero) ? null : new Parameter(cPtr, true);
    return ret;
  }

  public LocalParameter removeLocalParameter(long n) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_removeLocalParameter__SWIG_0(swigCPtr, n);
    LocalParameter ret = (cPtr == IntPtr.Zero) ? null : new LocalParameter(cPtr, true);
    return ret;
  }

  public Parameter removeParameter(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_removeParameter__SWIG_1(swigCPtr, sid);
    Parameter ret = (cPtr == IntPtr.Zero) ? null : new Parameter(cPtr, true);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public LocalParameter removeLocalParameter(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.KineticLaw_removeLocalParameter__SWIG_1(swigCPtr, sid);
    LocalParameter ret = (cPtr == IntPtr.Zero) ? null : new LocalParameter(cPtr, true);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.KineticLaw_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.KineticLaw_getElementName(swigCPtr);
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.KineticLaw_hasRequiredAttributes(swigCPtr);
    return ret;
  }

  public override bool hasRequiredElements() {
    bool ret = libsbmlPINVOKE.KineticLaw_hasRequiredElements(swigCPtr);
    return ret;
  }

}

}
