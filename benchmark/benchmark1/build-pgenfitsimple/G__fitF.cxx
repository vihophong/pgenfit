// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__fitF

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_fitF(void *p = 0);
   static void *newArray_fitF(Long_t size, void *p);
   static void delete_fitF(void *p);
   static void deleteArray_fitF(void *p);
   static void destruct_fitF(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::fitF*)
   {
      ::fitF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::fitF >(0);
      static ::ROOT::TGenericClassInfo 
         instance("fitF", ::fitF::Class_Version(), "include/fitF.hh", 29,
                  typeid(::fitF), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::fitF::Dictionary, isa_proxy, 4,
                  sizeof(::fitF) );
      instance.SetNew(&new_fitF);
      instance.SetNewArray(&newArray_fitF);
      instance.SetDelete(&delete_fitF);
      instance.SetDeleteArray(&deleteArray_fitF);
      instance.SetDestructor(&destruct_fitF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::fitF*)
   {
      return GenerateInitInstanceLocal((::fitF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::fitF*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_fitFbkg(void *p = 0);
   static void *newArray_fitFbkg(Long_t size, void *p);
   static void delete_fitFbkg(void *p);
   static void deleteArray_fitFbkg(void *p);
   static void destruct_fitFbkg(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::fitFbkg*)
   {
      ::fitFbkg *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::fitFbkg >(0);
      static ::ROOT::TGenericClassInfo 
         instance("fitFbkg", ::fitFbkg::Class_Version(), "include/fitF.hh", 93,
                  typeid(::fitFbkg), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::fitFbkg::Dictionary, isa_proxy, 4,
                  sizeof(::fitFbkg) );
      instance.SetNew(&new_fitFbkg);
      instance.SetNewArray(&newArray_fitFbkg);
      instance.SetDelete(&delete_fitFbkg);
      instance.SetDeleteArray(&deleteArray_fitFbkg);
      instance.SetDestructor(&destruct_fitFbkg);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::fitFbkg*)
   {
      return GenerateInitInstanceLocal((::fitFbkg*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::fitFbkg*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_fitFbkgflat(void *p = 0);
   static void *newArray_fitFbkgflat(Long_t size, void *p);
   static void delete_fitFbkgflat(void *p);
   static void deleteArray_fitFbkgflat(void *p);
   static void destruct_fitFbkgflat(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::fitFbkgflat*)
   {
      ::fitFbkgflat *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::fitFbkgflat >(0);
      static ::ROOT::TGenericClassInfo 
         instance("fitFbkgflat", ::fitFbkgflat::Class_Version(), "include/fitF.hh", 123,
                  typeid(::fitFbkgflat), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::fitFbkgflat::Dictionary, isa_proxy, 4,
                  sizeof(::fitFbkgflat) );
      instance.SetNew(&new_fitFbkgflat);
      instance.SetNewArray(&newArray_fitFbkgflat);
      instance.SetDelete(&delete_fitFbkgflat);
      instance.SetDeleteArray(&deleteArray_fitFbkgflat);
      instance.SetDestructor(&destruct_fitFbkgflat);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::fitFbkgflat*)
   {
      return GenerateInitInstanceLocal((::fitFbkgflat*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::fitFbkgflat*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr fitF::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *fitF::Class_Name()
{
   return "fitF";
}

//______________________________________________________________________________
const char *fitF::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitF*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int fitF::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitF*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *fitF::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitF*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *fitF::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitF*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr fitFbkg::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *fitFbkg::Class_Name()
{
   return "fitFbkg";
}

//______________________________________________________________________________
const char *fitFbkg::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitFbkg*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int fitFbkg::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitFbkg*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *fitFbkg::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitFbkg*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *fitFbkg::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitFbkg*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr fitFbkgflat::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *fitFbkgflat::Class_Name()
{
   return "fitFbkgflat";
}

//______________________________________________________________________________
const char *fitFbkgflat::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitFbkgflat*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int fitFbkgflat::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::fitFbkgflat*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *fitFbkgflat::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitFbkgflat*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *fitFbkgflat::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::fitFbkgflat*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void fitF::Streamer(TBuffer &R__b)
{
   // Stream an object of class fitF.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(fitF::Class(),this);
   } else {
      R__b.WriteClassBuffer(fitF::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_fitF(void *p) {
      return  p ? new(p) ::fitF : new ::fitF;
   }
   static void *newArray_fitF(Long_t nElements, void *p) {
      return p ? new(p) ::fitF[nElements] : new ::fitF[nElements];
   }
   // Wrapper around operator delete
   static void delete_fitF(void *p) {
      delete ((::fitF*)p);
   }
   static void deleteArray_fitF(void *p) {
      delete [] ((::fitF*)p);
   }
   static void destruct_fitF(void *p) {
      typedef ::fitF current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::fitF

//______________________________________________________________________________
void fitFbkg::Streamer(TBuffer &R__b)
{
   // Stream an object of class fitFbkg.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(fitFbkg::Class(),this);
   } else {
      R__b.WriteClassBuffer(fitFbkg::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_fitFbkg(void *p) {
      return  p ? new(p) ::fitFbkg : new ::fitFbkg;
   }
   static void *newArray_fitFbkg(Long_t nElements, void *p) {
      return p ? new(p) ::fitFbkg[nElements] : new ::fitFbkg[nElements];
   }
   // Wrapper around operator delete
   static void delete_fitFbkg(void *p) {
      delete ((::fitFbkg*)p);
   }
   static void deleteArray_fitFbkg(void *p) {
      delete [] ((::fitFbkg*)p);
   }
   static void destruct_fitFbkg(void *p) {
      typedef ::fitFbkg current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::fitFbkg

//______________________________________________________________________________
void fitFbkgflat::Streamer(TBuffer &R__b)
{
   // Stream an object of class fitFbkgflat.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(fitFbkgflat::Class(),this);
   } else {
      R__b.WriteClassBuffer(fitFbkgflat::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_fitFbkgflat(void *p) {
      return  p ? new(p) ::fitFbkgflat : new ::fitFbkgflat;
   }
   static void *newArray_fitFbkgflat(Long_t nElements, void *p) {
      return p ? new(p) ::fitFbkgflat[nElements] : new ::fitFbkgflat[nElements];
   }
   // Wrapper around operator delete
   static void delete_fitFbkgflat(void *p) {
      delete ((::fitFbkgflat*)p);
   }
   static void deleteArray_fitFbkgflat(void *p) {
      delete [] ((::fitFbkgflat*)p);
   }
   static void destruct_fitFbkgflat(void *p) {
      typedef ::fitFbkgflat current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::fitFbkgflat

namespace {
  void TriggerDictionaryInitialization_libfitF_Impl() {
    static const char* headers[] = {
"/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh",
0
    };
    static const char* includePaths[] = {
"/opt/cernroot/root_v6.08.00/include",
"/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple",
"/opt/cernroot/root_v6.08.00/include",
"/data03/users/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/build-pgenfitsimple/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libfitF dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh")))  fitF;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh")))  fitFbkg;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh")))  fitFbkgflat;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libfitF dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/phong/projects/brikenmacros/fitting/pgenfit/pgenfitv5/benchmark/benchmark1/pgenfitsimple/include/fitF.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"fitF", payloadCode, "@",
"fitFbkg", payloadCode, "@",
"fitFbkgflat", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libfitF",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libfitF_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libfitF_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libfitF() {
  TriggerDictionaryInitialization_libfitF_Impl();
}
