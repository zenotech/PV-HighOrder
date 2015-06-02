//Paraview plugin for the visualization of high order fields
//Sebastien Blaise, 2012

#ifndef _vtkGalerkin_h
#define _vtkGalerkin_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkEdgeTable.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <map>
#include <vector>
#include <iostream>
#include <stdlib.h>



class vtkHighOrder : public vtkUnstructuredGridAlgorithm
{
  
  class shapeFunctions;
  class recurShape;
  
  
  //Useful data for splitting and writing elements
  class adaptData
  {
  public:
    int isTwoLevel;
    int iError;
    float sqAvg;
    float tolerance;
    vtkFloatArray **arrays;
    int *dim; 
    int nbFields;
    vtkUnstructuredGrid *out;
    vtkPoints* newPts;
    vtkHighOrder* plugin;
    adaptData(int nbPointArrays, int* dim_, std::string *fieldsNames);
    ~adaptData();
  };

  class adaptDataElem
  {
  public:
    float ***dof;
    float **dofCoord;
    int nbShapeFct;
  };
  
  //A general entity to be recursively split
  class adapt_entity
  {
  protected:
    void writeEntity(float **valShape,float **valShapeLin);
    //    virtual void recursion_split(int leftlevels=0)=0;
    int elem_type;
    int npoints,nsubEntities;
    float **points;
    adaptData *data;
    adaptDataElem *dataElem;
    void getAvg(recurShape *tree, float *avg);
  public:
    adapt_entity(adaptData *data_,adaptDataElem *dataElem_,int npoints_,int nsubEntities_,int elem_type_,
		 recurShape *tree, recurShape *treeLin, int leftlevels);
    ~adapt_entity();
  };
    
 private:
  //Genreal tools
  vtkUnstructuredGrid *in;
  vtkUnstructuredGrid *out;
  
  int levelMax;
  std::string errorField;
  float relTolerance;
  float isTwoLevel;

  int nbElems;
  vtkPoints* newPts;
  
  recurShape *getQuadTree(int nbpoints, int maxlevel);
  recurShape *getTriangleTree(int nbpoints, int maxlevel);
  recurShape *getTetrahedronTree(int nbpoints, int maxlevel);
  recurShape *getHexahedronTree(int nbpoints, int maxlevel);
    
 protected:
  
  void printErrorMessage(std::string msg){
    vtkErrorMacro(<< msg);
  }
  
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int, vtkInformation *info);
  void subdivideAll();
  
  vtkTypeMacro(vtkHighOrder,vtkUnstructuredGridAlgorithm);
  vtkHighOrder(){};
  ~vtkHighOrder(){};
  
 public:
  //Interface
  void SetLevelMax(int levelMax);
  void SetErrorField(int a, int b, int c, int d, const char* errorField);
  void SetRelTolerance(float relTolerance_);
  void SetTwoLevelErr(int isTwoLevel);

  static vtkHighOrder *New();
  
  
};


#endif
