// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libTAUn_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#include "taun_vars.h"
#include "TauData.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TauData_Dictionary();
   static void TauData_TClassManip(TClass*);
   static void *new_TauData(void *p = 0);
   static void *newArray_TauData(Long_t size, void *p);
   static void delete_TauData(void *p);
   static void deleteArray_TauData(void *p);
   static void destruct_TauData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TauData*)
   {
      ::TauData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TauData));
      static ::ROOT::TGenericClassInfo 
         instance("TauData", "TauData.h", 9,
                  typeid(::TauData), DefineBehavior(ptr, ptr),
                  &TauData_Dictionary, isa_proxy, 0,
                  sizeof(::TauData) );
      instance.SetNew(&new_TauData);
      instance.SetNewArray(&newArray_TauData);
      instance.SetDelete(&delete_TauData);
      instance.SetDeleteArray(&deleteArray_TauData);
      instance.SetDestructor(&destruct_TauData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TauData*)
   {
      return GenerateInitInstanceLocal((::TauData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TauData*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TauData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TauData*)0x0)->GetClass();
      TauData_TClassManip(theClass);
   return theClass;
   }

   static void TauData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_TauData(void *p) {
      return  p ? new(p) ::TauData : new ::TauData;
   }
   static void *newArray_TauData(Long_t nElements, void *p) {
      return p ? new(p) ::TauData[nElements] : new ::TauData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TauData(void *p) {
      delete ((::TauData*)p);
   }
   static void deleteArray_TauData(void *p) {
      delete [] ((::TauData*)p);
   }
   static void destruct_TauData(void *p) {
      typedef ::TauData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TauData

namespace {
  void TriggerDictionaryInitialization_libTAUn_Dict_Impl() {
    static const char* headers[] = {
"taun_vars.h",
"TauData.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include/root",
"/home/grammer/anal/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TauData.h")))  TauData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "taun_vars.h"
#include "TauData.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TauData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libTAUn_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libTAUn_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libTAUn_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libTAUn_Dict() {
  TriggerDictionaryInitialization_libTAUn_Dict_Impl();
}
