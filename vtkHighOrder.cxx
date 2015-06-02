//Paraview plugin for the visualization of high order fields
//Sebastien Blaise, 2012


#include "vtkHighOrder.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolygon.h"
#include "vtkStreamingDemandDrivenPipeline.h"


#include <sstream>
#include <set>
#include <vector>

vtkStandardNewMacro(vtkHighOrder);

#include "shapeFunctions.h"

int vtkHighOrder::RequestData(
			      vtkInformation *vtkNotUsed(request),
			      vtkInformationVector **inputVector,
			      vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  in=vtkUnstructuredGrid::SafeDownCast(this->GetInput());
  nbElems = in->GetNumberOfCells();
  
  out=this->GetOutput();
  out->Initialize();
  out->Allocate();
  
  //Output points 
  newPts = vtkPoints::New();
  
  subdivideAll();
  
  out->SetPoints(newPts);
  newPts->Delete();
  return 1;
}

int vtkHighOrder::FillInputPortInformation(int, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
}

void vtkHighOrder::subdivideAll(){
  
  vtkIdType nbPointArrays = in->GetPointData()->GetNumberOfArrays();
  vtkIdType nbCellArrays = in->GetCellData()->GetNumberOfArrays();

  if(nbPointArrays == 0)
  {
	  vtkErrorMacro(<< "No PointData arrays found");
  }
  if(nbCellArrays == 0)
  {
	  vtkErrorMacro(<< "No CellData arrays found");
  }

  // Create a set of all variables with HO solution pts
  std::set<std::string> HOvariables;
  for(vtkIdType j=0; j < nbCellArrays; ++j)
  {
	  std::string name = in->GetCellData()->GetArrayName(j);
	  std::size_t pos;
	  if((pos = name.find("_HOsol_")) !=std::string::npos)
	  {
		  name = name.substr(0,pos);
		  // Check if point data array exists
		  if(vtkFloatArray::SafeDownCast(in->GetPointData()->GetArray(name.c_str())))
			  HOvariables.insert(name);
	  }
  }

  // Check
  int nbAddDof=0;
  for(vtkIdType j = 0; j < nbPointArrays; j++)
    {
      std::string namePt=in->GetPointData()->GetArrayName(j);

      if(HOvariables.find(namePt) != HOvariables.end())
      {
		  int iAddDof=-1;
		  while (true){
			iAddDof++;
			std::stringstream nameCSol;
			nameCSol << namePt << "_HOsol_" << iAddDof;
			if(!vtkFloatArray::SafeDownCast(in->GetCellData()->GetArray(nameCSol.str().c_str())))
			  break;
		  }
		  if (nbAddDof!=0 && iAddDof!=nbAddDof)
			  vtkErrorMacro(<< "All PointData fields should have the same polynomial order");
		  nbAddDof=iAddDof;
      }
    }
  
  std::map<std::string,std::vector<vtkFloatArray*> > arraysDof;
  std::vector<vtkFloatArray*> arraysCoord;

  //vtkFloatArray *arraysDof[nbPointArrays][nbAddDof];
  //vtkFloatArray *arraysCoord[nbAddDof];
  nbPointArrays = 0;
  for(std::set<std::string>::iterator itr=HOvariables.begin();itr!=HOvariables.end();++itr)
  {
	  for (int i=0; i<nbAddDof; i++)
	  {
	  std::stringstream nameCSol;
	  nameCSol << *itr << "_HOsol_" << i;
	  arraysDof[*itr].push_back(
			  vtkFloatArray::SafeDownCast(in->GetCellData()->GetArray(nameCSol.str().c_str())));
	  }
	  nbPointArrays++;
  }

  bool highOrderGeo=true;
  for (int i=0; i<nbAddDof; i++)
    {
      std::stringstream nameCCoord;
      nameCCoord << "HOcoord_" << i;
      arraysCoord.push_back(vtkFloatArray::SafeDownCast(in->GetCellData()->GetArray(nameCCoord.str().c_str())));
      if (!arraysCoord[i]) {
    	  highOrderGeo=false;
      }	  
    }
  
  std::vector<std::string> fieldsNames(nbPointArrays);
  int dim[nbPointArrays];
  nbPointArrays=0;
  for(std::set<std::string>::iterator itr=HOvariables.begin();itr!=HOvariables.end();++itr)
  {
  //for(vtkIdType j = 0; j < nbPointArrays; j++){
    dim[nbPointArrays]=arraysDof[*itr][0]->GetNumberOfComponents();
    fieldsNames[nbPointArrays]=*itr;
    nbPointArrays++;
  }
  adaptData data(nbPointArrays,dim,fieldsNames);
  data.out=out;
  data.newPts=newPts;
  data.plugin=this;
  data.tolerance=relTolerance;
  data.isTwoLevel=isTwoLevel;

  for(vtkIdType j = 0; j < nbPointArrays; j++){
    if (fieldsNames[j]==errorField)
      data.iError=j;
  }

  //Compute the average value of errorField
  vtkFloatArray *errArray=vtkFloatArray::SafeDownCast(in->GetPointData()->GetArray(fieldsNames[data.iError].c_str()));
  float val[dim[data.iError]];

  float sqAvg=0;
  for (int i=0; i<in->GetNumberOfPoints(); i++){
    errArray->GetTupleValue(i,val);
    for (int j=0; j<dim[data.iError]; j++){
      sqAvg+=val[j]*val[j];
    }
  }
  float val2[dim[data.iError]];
  for (int k=0; k<nbAddDof; k++){
    for (int i=0; i<in->GetNumberOfCells(); i++){

      arraysDof[fieldsNames[data.iError]][k]->GetTupleValue(i,val2);

      for (int j=0; j<dim[data.iError]; j++){
    	  sqAvg+=val2[j]*val2[j];
      }
    }
  }
  sqAvg /= (in->GetNumberOfPoints()+in->GetNumberOfCells()*nbAddDof);
  data.sqAvg=sqAvg;

  int iElem;

  //We build the recurrence trees (serial because singleton)
  for (iElem=0; iElem<nbElems; iElem++)
    {
      vtkIdList *idsIn = vtkIdList::New();
      in->GetCellPoints(iElem,idsIn);
      
      int nbNodesByCell=idsIn->GetNumberOfIds();
      int cellType=in->GetCellType(iElem);

      switch (cellType){
      case VTK_TRIANGLE:
		getTriangleTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
		if (!highOrderGeo)
		  getTriangleTree(nbNodesByCell, levelMax);
		break;
      case VTK_QUAD:
        getQuadTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
        if (!highOrderGeo)
          getQuadTree(nbNodesByCell, levelMax);
        break;
      case VTK_TETRA:
        getTetrahedronTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
        if (!highOrderGeo)
          getTetrahedronTree(nbNodesByCell, levelMax);
        break;
      case VTK_HEXAHEDRON:
		getHexahedronTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
		if (!highOrderGeo)
		  getHexahedronTree(nbNodesByCell, levelMax);
		break;
      default:
	vtkErrorMacro(<< "Element type unknown ("<<nbNodesByCell<<") nodes"); 
      }
    }

#pragma omp parallel for
  for (iElem=0; iElem<nbElems; iElem++)
    {
      
      vtkIdType np;
      vtkIdList *idsIn = vtkIdList::New();;
      in->GetCellPoints(iElem,idsIn);
      
      int nbNodesByCell=idsIn->GetNumberOfIds();
      int cellType=in->GetCellType(iElem);
      
      float **dof[nbPointArrays];
      int sizeDofCoord=highOrderGeo?(nbNodesByCell+nbAddDof):nbNodesByCell;
      float *dofCoord[sizeDofCoord];
      for (int i=0; i<sizeDofCoord; i++)
    	  dofCoord[i]=new float[3];
      for (int i=0; i<nbPointArrays; i++)
    	  dof[i]=new float*[nbNodesByCell+nbAddDof];
      
      //P1 data
      for (int iNode=0; iNode<nbNodesByCell; iNode++){
    	double v[3];
		in->GetPoint(idsIn->GetId(iNode),v);
		for (int i=0; i<3; i++)
			dofCoord[iNode][i] = v[i];
		for (int j=0; j<nbPointArrays; j++){
		  vtkFloatArray *ptarray=vtkFloatArray::SafeDownCast(in->GetPointData()->GetArray(fieldsNames[j].c_str()));
		  dof[j][iNode]=new float[dim[j]];
		  ptarray->GetTupleValue(idsIn->GetId(iNode),dof[j][iNode]);
		}
      }
      //High Order extension
      for (int i=0; i<nbAddDof; i++){
	for (int j=0; j<nbPointArrays; j++){
	  dof[j][i+nbNodesByCell]=new float[dim[j]];
	  arraysDof[fieldsNames[j]][i]->GetTupleValue (iElem, dof[j][i+nbNodesByCell]);
	}
	if (highOrderGeo)
	  arraysCoord[i]->GetTupleValue (iElem, dofCoord[i+nbNodesByCell]);
      }

      adaptDataElem dataElem;
      dataElem.dof=dof;
      dataElem.dofCoord=dofCoord;
      dataElem.nbShapeFct=nbNodesByCell+nbAddDof;
      
      recurShape *tree,*treeLin;
      switch (cellType){
      case VTK_TRIANGLE:
	{
	  //#pragma omp critical(buildTreeTRI)
	  {
	    tree=getTriangleTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
	    if (!highOrderGeo)
	      treeLin=getTriangleTree(nbNodesByCell, levelMax);
	    else
	      treeLin=NULL;
	  }
	  if (tree)
	    adapt_entity triangle(&data,&dataElem,3,4,VTK_TRIANGLE,tree,treeLin,levelMax);
	  break;
	}
      case VTK_QUAD:
	{
          //#pragma omp critical(buildTreeQUAD)
          {
            tree=getQuadTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
            if (!highOrderGeo)
              treeLin=getQuadTree(nbNodesByCell, levelMax);
            else
              treeLin=NULL;
          }
          if (tree)
            adapt_entity quad(&data,&dataElem,4,4,VTK_QUAD,tree,treeLin,levelMax);
          break;
        }
      case VTK_TETRA:
        {
          //#pragma omp critical(buildTreeTET)
          {
            tree=getTetrahedronTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
            if (!highOrderGeo)
              treeLin=getTetrahedronTree(nbNodesByCell, levelMax);
            else
              treeLin=NULL;
          }
          if (tree)
            adapt_entity tetrahedron(&data,&dataElem,4,8,VTK_TETRA,tree,treeLin,levelMax);
          break;
        }
      case VTK_HEXAHEDRON:
	{
	  //#pragma omp critical(buildTreeHEX)
	  {
	    tree=getHexahedronTree(nbNodesByCell+nbAddDof, levelMax+((isTwoLevel==1)?1:0));
	    if (!highOrderGeo)
	      treeLin=getHexahedronTree(nbNodesByCell, levelMax);
	    else
	      treeLin=NULL;
	  }
	  if (tree)
	    adapt_entity hexahedron(&data,&dataElem,8,8,VTK_HEXAHEDRON,tree,treeLin,levelMax);
	  break;
	}
      default:
	vtkErrorMacro(<< "Element type unknown ("<<cellType<<") nodes"); 
      }
      
      
      //Clearing memory
      for (int i=0; i<sizeDofCoord; i++)
	delete[] dofCoord[i];
      for (int j=0; j<nbPointArrays; j++){
	for (int i=0; i<nbNodesByCell+nbAddDof; i++){
	  delete[] dof[j][i];
	}
	delete[] dof[j];
      }
      idsIn->Delete();
    }

}

void vtkHighOrder::SetLevelMax(int levelMax_){
  levelMax=levelMax_;
  this->Modified();
}

void vtkHighOrder::SetErrorField(int a, int b, int c, int d, const char* errorField_){
  errorField=errorField_;
  this->Modified();
}

void vtkHighOrder::SetRelTolerance(float relTolerance_){
  relTolerance=relTolerance_;
  this->Modified();
}

void vtkHighOrder::SetTwoLevelErr(int isTwoLevel_){
  isTwoLevel=isTwoLevel_;
  this->Modified();
}

vtkHighOrder::adaptData::adaptData(int nbPointArrays, int* dim_, std::vector<std::string> &fieldsNames){
  nbFields=nbPointArrays;
  dim=dim_;
  arrays = new vtkFloatArray*[nbFields];
  for (int i=0; i<nbFields; i++){
    arrays[i] = vtkFloatArray::New();
    arrays[i]->SetNumberOfComponents(dim[i]);
    arrays[i]->SetName(fieldsNames[i].c_str());
  }
}
vtkHighOrder::adaptData::~adaptData(){
  //Pass the data arrays to output vtkUnstructuredGrid before deleting 
  for (int i=0; i<nbFields; i++){
      out->GetPointData()->AddArray(arrays[i]);
      arrays[i]->Delete();
  }
  delete[] arrays;
}


vtkHighOrder::adapt_entity::adapt_entity(adaptData *data_,adaptDataElem *dataElem_, int npoints_,int nsubEntities_,int elem_type_,
					 recurShape *tree, recurShape *treeLin, int leftlevels)
  :data(data_),dataElem(dataElem_),npoints(npoints_),nsubEntities(nsubEntities_),elem_type(elem_type_)
{
  // 1 - check the error to know if children are needed
  bool needRefine=false;
  if (data->tolerance<0)
    needRefine=true;
  else if (tree->getChild(0)){
    float avgErr[data->dim[data->iError]];
    float avg[data->dim[data->iError]];
    //Two-levels check
    if (tree->getChild(0)->getChild(0) && data->isTwoLevel==1){
      float sqerror=0;
      for (int i=0; i<nsubEntities; i++){
	getAvg(tree->getChild(i),avgErr); //Average on a child
	for (int j=0; j<nsubEntities; j++){
	  getAvg(tree->getChild(i)->getChild(j),avg); //Average on a sub-child
	  for (int k=0; k<data->dim[data->iError]; k++){
	    avgErr[k]-=avg[k]/nsubEntities;
	  }
	}
	for (int j=0; j<data->dim[data->iError]; j++)
	  sqerror+=avgErr[j]*avgErr[j];
	if (sqerror>data->tolerance*data->tolerance*data->sqAvg)
	  needRefine=true;
      }
    }
    //One-level check
    float sqerror=0;
    getAvg(tree,avgErr); //Average on current element
    for (int i=0; i<nsubEntities; i++){
      getAvg(tree->getChild(i),avg); //Average on a child
      for (int j=0; j<data->dim[data->iError]; j++){
	avgErr[j]-=avg[j]/nsubEntities;
      }
    }
    for (int j=0; j<data->dim[data->iError]; j++)
      sqerror+=avgErr[j]*avgErr[j];
    if (sqerror>data->tolerance*data->tolerance*data->sqAvg)
      needRefine=true;
  }
  
  // 2a - no children
  if (!needRefine || leftlevels <= 0 || !tree->getChild(0)){
    if (leftlevels > 0 && !tree->getChild(0)){
      data->plugin->printErrorMessage("Maximum level of available refinements reached");
    }
    double **valShape=tree->getValShape();
    double **valShapeLin=NULL;
    if (treeLin)
      valShapeLin=treeLin->getValShape();
    //#pragma omp critical(writeEntity)
    {
    writeEntity(valShape,valShapeLin);
    }
    return;
  }
  
  // 2b - recursively build children
  for (int i=0; i<nsubEntities; i++){
    adapt_entity(data,dataElem,npoints_,nsubEntities_,elem_type_,tree->getChild(i),(treeLin?treeLin->getChild(i):NULL),leftlevels-1);
  }
  
}

vtkHighOrder::adapt_entity::~adapt_entity(){
}

void vtkHighOrder::adapt_entity::writeEntity(double **valShape,double **valShapeLin){
  float points[npoints][3];
  //Coordinates
  vtkIdList *ids = vtkIdList::New();

  for (int k=0; k<npoints; k++){
    for (int j=0; j<3; j++){
      points[k][j]=0;
      if (!valShapeLin){
	for (int i=0; i<dataElem->nbShapeFct; i++){
	  points[k][j]=points[k][j]+valShape[k][i]*dataElem->dofCoord[i][j];
	}
      }else{
	for (int i=0; i<npoints; i++){
	  points[k][j]=points[k][j]+valShapeLin[k][i]*dataElem->dofCoord[i][j];
	}
      }
    }
  }
  int maxdim=0;
  for (int iF=0; iF<data->nbFields; iF++)
    maxdim=std::max(maxdim,data->dim[iF]);
  float val[npoints][data->nbFields][maxdim];
  //Data
  for (int k=0; k<npoints; k++){
    for (int iF=0; iF<data->nbFields; iF++){
      for (int j=0; j<data->dim[iF]; j++){
	val[k][iF][j]=0;
	for (int i=0; i<dataElem->nbShapeFct; i++){
	  val[k][iF][j]+=valShape[k][i]*dataElem->dof[iF][i][j];
	}
      }
    }
  }
  #pragma omp critical(writeElemVtk)
  {
    for (int k=0; k<npoints; k++){
      for (int iF=0; iF<data->nbFields; iF++){
	data->arrays[iF]->InsertNextTupleValue(val[k][iF]);
      }
      ids->InsertNextId(data->newPts->InsertNextPoint(points[k]));
    }
    data->out->InsertNextCell(elem_type,ids);   
  }
  ids->Delete();
}


void vtkHighOrder::adapt_entity::getAvg(recurShape *tree, float *avg){
  double **valShape=tree->getValShape();
  for (int j=0; j<data->dim[data->iError]; j++){
    avg[j]=0;
    for (int k=0; k<npoints; k++){
      for (int i=0; i<dataElem->nbShapeFct; i++){
	avg[j]+=valShape[k][i]*dataElem->dof[data->iError][i][j];
      }
    }
    avg[j]=avg[j]/npoints;
  }
}


