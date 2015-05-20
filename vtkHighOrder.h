//Paraview plugin for the visualization of high order fields
//Sebastien Blaise, 2012

#ifndef _vtkGalerkin_h
#define _vtkGalerkin_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkEdgeTable.h"
#include "vtkDoubleArray.h"
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
    double sqAvg;
    double tolerance;
    vtkDoubleArray **arrays;
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
    double ***dof;
    double **dofCoord; 
    int nbShapeFct;
  };
  
  //A general entity to be recursively split
  class adapt_entity
  {
  protected:
    void writeEntity(double **valShape,double **valShapeLin);
    //    virtual void recursion_split(int leftlevels=0)=0;
    int elem_type;
    int npoints,nsubEntities;
    double **points;
    adaptData *data;
    adaptDataElem *dataElem;
    void getAvg(recurShape *tree, double *avg);
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
  double relTolerance;
  double isTwoLevel;

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
  void SetRelTolerance(double relTolerance_);
  void SetTwoLevelErr(int isTwoLevel);

  static vtkHighOrder *New();
  
  
};


#endif
