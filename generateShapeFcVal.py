from gmshpy import *
from math import *
from array import *
import StringIO
import os, sys

#QUAD
typesQ=[MSH_QUA_1,MSH_QUA_4,MSH_QUA_9,MSH_QUA_16,MSH_QUA_25,MSH_QUA_36,MSH_QUA_49,MSH_QUA_64,MSH_QUA_81,MSH_QUA_100,MSH_QUA_121]
nQ=[1,4,9,16,25,36,49,64,81,100,121]

#TRI
typesT=[MSH_TRI_1,MSH_TRI_3,MSH_TRI_6,MSH_TRI_10,MSH_TRI_15,MSH_TRI_21,MSH_TRI_28,MSH_TRI_36,MSH_TRI_45,MSH_TRI_55,MSH_TRI_66]
nT=[1,3,6,10,15,21,28,36,45,55,66]

#HEX
typesH=[MSH_HEX_8,MSH_HEX_27,MSH_HEX_64,MSH_HEX_125,MSH_HEX_216,MSH_HEX_343,MSH_HEX_512,MSH_HEX_729,MSH_HEX_1000]
nH=[8,27,64,125,216,343,512,729,1000]

#ALL
geomTypes=["Quad","Triangle","Hexahedron"]
allGmshTypes=[typesQ,typesT,typesH]
allN=[nQ,nT,nH]
nbChilds=[4,4,8]

def Entity(elemType,iGeom,fsLin,points,pointsLin,pointsIso,children):
    for i in range(nbChilds[iGeom]):
        for j in range(points.size1()):
            for k in range(points.size2()):
                pointsIso[i].set(j,k,0)

        shapeLin=fullMatrixDouble(points.size1(),pointsLin.size1());
        fsLin.f(points,shapeLin)
        for iPts in range(points.size1()):
            for iDim in range(points.size2()):
                for iShapeLin in range(pointsLin.size1()):
                    pointsIso[i].set(iPts,iDim,pointsIso[i](iPts,iDim)+shapeLin(iPts,iShapeLin)*children[i][iShapeLin](0,iDim))

    print >>classdef, ("   static const double childVal[%i][%i][%i];")%(nbChilds[iGeom],points.size1(),points.size1())
    print >>outclass, ("const double vtkHighOrder::shapeFunctions::%s%i::childVal[%i][%i][%i]={")%(geomTypes[iGeom],elemType,nbChilds[iGeom],points.size1(),points.size1()),
    for i in range(len(pointsIso)):
        shapeF=fullMatrixDouble(points.size1(),points.size1())
        fs.f(pointsIso[i],shapeF)
        for iPts in range(points.size1()):
            for iShape in range(points.size1()):
                print >>outclass, shapeF(iPts,iShape),
                if (iPts+iShape+i<len(pointsIso)-1+2*(points.size1()-1)):
                    print >>outclass, ",",
    print >>outclass, "};"

def Quad(elemType):
    fsLin = BasisFactory.create(MSH_QUA_4)
    pointsLin=fsLin.points

    points=fs.points
    p1=fullMatrixDouble(1,points.size2())
    p2=fullMatrixDouble(1,points.size2())
    p3=fullMatrixDouble(1,points.size2())
    p4=fullMatrixDouble(1,points.size2())
    for j in range (pointsLin.size2()):
        p1.set(0,j,pointsLin(0,j));
        p2.set(0,j,pointsLin(1,j));
        p3.set(0,j,pointsLin(2,j));
        p4.set(0,j,pointsLin(3,j));
    p12=fullMatrixDouble(1,points.size2());
    p23=fullMatrixDouble(1,points.size2());
    p34=fullMatrixDouble(1,points.size2());
    p14=fullMatrixDouble(1,points.size2());
    pc=fullMatrixDouble(1,points.size2());
    for i in range(points.size2()):
        p12.set(0,i,(p1(0,i) + p2(0,i))*0.5)
        p23.set(0,i,(p2(0,i) + p3(0,i))*0.5)
        p34.set(0,i,(p3(0,i) + p4(0,i))*0.5)
        p14.set(0,i,(p1(0,i) + p4(0,i))*0.5)
        pc.set(0,i,(p1(0,i) + p2(0,i) + p3(0,i) + p4(0,i))*0.25)

    child1=[p1,p12,pc,p14]
    child2=[p2,p23,pc,p12]
    child3=[p3,p34,pc,p23]
    child4=[p4,p14,pc,p34]
    children=[child1,child2,child3,child4]

    pointsIso=[fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2())]
    Entity(elemType,0,fsLin,points,pointsLin,pointsIso,children)

def Triangle(elemType):
    points=fs.points
    fsLin = BasisFactory.create(MSH_TRI_3)
    pointsLin=fsLin.points

    p1=fullMatrixDouble(1,points.size2())
    p2=fullMatrixDouble(1,points.size2())
    p3=fullMatrixDouble(1,points.size2())
    for j in range (pointsLin.size2()):
        p1.set(0,j,pointsLin(0,j));
        p2.set(0,j,pointsLin(1,j));
        p3.set(0,j,pointsLin(2,j));
    p12=fullMatrixDouble(1,points.size2());
    p13=fullMatrixDouble(1,points.size2());
    p23=fullMatrixDouble(1,points.size2());
    for i in range(points.size2()):
        p12.set(0,i,(p1(0,i) + p2(0,i))*0.5)
        p13.set(0,i,(p1(0,i) + p3(0,i))*0.5)
        p23.set(0,i,(p2(0,i) + p3(0,i))*0.5)

    child1=[p1,p12,p13]
    child2=[p2,p23,p12]
    child3=[p3,p13,p23]
    child4=[p12,p23,p13]
    children=[child1,child2,child3,child4]

    pointsIso=[fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2())]
    Entity(elemType,1,fsLin,points,pointsLin,pointsIso,children)

def Hexahedron(elemType):
    fsLin = BasisFactory.create(MSH_HEX_8)
    pointsLin=fsLin.points

    points=fs.points
    p0=fullMatrixDouble(1,points.size2())
    p1=fullMatrixDouble(1,points.size2())
    p2=fullMatrixDouble(1,points.size2())
    p3=fullMatrixDouble(1,points.size2())
    p4=fullMatrixDouble(1,points.size2())
    p5=fullMatrixDouble(1,points.size2())
    p6=fullMatrixDouble(1,points.size2())
    p7=fullMatrixDouble(1,points.size2())
    for j in range (pointsLin.size2()):
        p0.set(0,j,pointsLin(0,j))
        p1.set(0,j,pointsLin(1,j))
        p2.set(0,j,pointsLin(2,j))
        p3.set(0,j,pointsLin(3,j))
        p4.set(0,j,pointsLin(4,j))
        p5.set(0,j,pointsLin(5,j))
        p6.set(0,j,pointsLin(6,j))
        p7.set(0,j,pointsLin(7,j))
        p01=fullMatrixDouble(1,points.size2())
        p12=fullMatrixDouble(1,points.size2())
        p23=fullMatrixDouble(1,points.size2())
        p03=fullMatrixDouble(1,points.size2())
        p45=fullMatrixDouble(1,points.size2())
        p56=fullMatrixDouble(1,points.size2())
        p67=fullMatrixDouble(1,points.size2())
        p47=fullMatrixDouble(1,points.size2())
        p04=fullMatrixDouble(1,points.size2())
        p15=fullMatrixDouble(1,points.size2())
        p26=fullMatrixDouble(1,points.size2())
        p37=fullMatrixDouble(1,points.size2())
        p0145=fullMatrixDouble(1,points.size2())
        p1256=fullMatrixDouble(1,points.size2())
        p2367=fullMatrixDouble(1,points.size2())
        p0347=fullMatrixDouble(1,points.size2())
        p4756=fullMatrixDouble(1,points.size2())
        p0312=fullMatrixDouble(1,points.size2())
        pc=fullMatrixDouble(1,points.size2())

    for i in range(points.size2()):
        p01.set(0,i,((p0(0,i) + p1(0,i)) * 0.5))
        p12.set(0,i,((p1(0,i) + p2(0,i)) * 0.5))
        p23.set(0,i,((p2(0,i) + p3(0,i)) * 0.5))
        p03.set(0,i,((p3(0,i) + p0(0,i)) * 0.5))
        p45.set(0,i,((p4(0,i) + p5(0,i)) * 0.5))
        p56.set(0,i,((p5(0,i) + p6(0,i)) * 0.5))
        p67.set(0,i,((p6(0,i) + p7(0,i)) * 0.5))
        p47.set(0,i,((p7(0,i) + p4(0,i)) * 0.5))
        p04.set(0,i,((p4(0,i) + p0(0,i)) * 0.5))
        p15.set(0,i,((p5(0,i) + p1(0,i)) * 0.5))
        p26.set(0,i,((p6(0,i) + p2(0,i)) * 0.5))
        p37.set(0,i,((p7(0,i) + p3(0,i)) * 0.5))
        p0145.set(0,i,((p45(0,i) + p01(0,i)) * 0.5))
        p1256.set(0,i,((p12(0,i) + p56(0,i)) * 0.5))
        p2367.set(0,i,((p23(0,i) + p67(0,i)) * 0.5))
        p0347.set(0,i,((p03(0,i) + p47(0,i)) * 0.5))
        p4756.set(0,i,((p47(0,i) + p56(0,i)) * 0.5))
        p0312.set(0,i,((p03(0,i) + p12(0,i)) * 0.5))
        pc.set(0,i,((p0(0,i) + p1(0,i) + p2(0,i) + p3(0,i) + p4(0,i) + p5(0,i) + p6(0,i) + p7(0,i)) * 0.125))


    child1=[p0, p01, p0312, p03, p04, p0145, pc, p0347]
    child2=[p01, p0145, p15, p1, p0312, pc, p1256, p12]
    child3=[p04, p4, p45, p0145, p0347, p47, p4756, pc]
    child4=[p0145, p45, p5, p15, pc, p4756, p56, p1256]
    child5=[p0347, p47, p4756, pc, p37, p7, p67, p2367]
    child6=[pc, p4756, p56, p1256, p2367, p67, p6, p26]
    child7=[p03, p0347, pc, p0312, p3, p37, p2367, p23]
    child8=[p0312, pc, p1256, p12, p23, p2367, p26, p2]
    children=[child1,child2,child3,child4,child5,child6,child7,child8]

    pointsIso=[fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2()),
               fullMatrixDouble(points.size1(),points.size2())]
    Entity(elemType,2,fsLin,points,pointsLin,pointsIso,children)


allFunctions=[Quad,Triangle,Hexahedron]

print "//Paraview plugin for the visualization of high order fields"
print "//Description of the shape functions"
print "//Sebastien Blaise, 2012"
print ""
print "//This code was generated automatically"
print "//Memory is not freed but it is not necessary (singleton class)"
print "#include <stdlib.h>"

print "class vtkHighOrder::shapeFunctions {"
print "public:"
for i in range(len(geomTypes)):
    for j in range(len(allGmshTypes[i])):
        print ("   class %s%i;")%(geomTypes[i],allN[i][j])
print "};"
print "class vtkHighOrder::recurShape {"
print "   double **valShape;"
print "   recurShape **children;"
print "public:"
print "   recurShape(double **valShape_, recurShape **children_){"
print "      valShape=valShape_;"
print "      children=children_;"
print "   }"
print "   void setData(double **valShape_, recurShape **children_){"
print "      valShape=valShape_;"
print "      children=children_;"
print "   }"
print "   double **getValShape(){"
print "      return valShape;"
print "   }"
print "   recurShape **getChildren(){"
print "      return children;"
print "   }"
print "   recurShape *getChild(int i){"
print "      return children[i];"
print "   }"
print "};"

p1=fullMatrixDouble(1,2)
p2=fullMatrixDouble(1,2)
p3=fullMatrixDouble(1,2)
p4=fullMatrixDouble(1,2)
p1.set(0,0,-1)
p1.set(0,1,-1)
p2.set(0,0,1)
p2.set(0,1,-1)
p3.set(0,0,1)
p3.set(0,1,1)
p4.set(0,0,-1)
p4.set(0,1,1)
classdef = StringIO.StringIO()
classconstr = StringIO.StringIO()
outclass = StringIO.StringIO()

for iGeom in range(len(geomTypes)):
    for i in range(len(allGmshTypes[iGeom])):
        classdef.truncate(0)
        classconstr.truncate(0)
        outclass.truncate(0)
        fs = BasisFactory.create(allGmshTypes[iGeom][i])
        allFunctions[iGeom](allN[iGeom][i])
        print ("class vtkHighOrder::shapeFunctions::%s%i {")%(geomTypes[iGeom],allN[iGeom][i])
        print ("private:");
        print ("   static int maxLevel;");
        print classdef.getvalue(),
        print ("   static recurShape *tree;")
        print ("   void  recurChild(int leftlevel, recurShape &parent){")
        print ("       if (leftlevel-- > 0){");
        print ("          for (int i=0; i<%i; i++){")%nbChilds[iGeom]
        print ("              if (parent.getChild(i)==NULL){ //They may be already built as we sometimes extend the level");
        print ("                    double **valShape=new double*[%i];")%(fs.points.size1())
        print ("                 for (int j=0; j<%i; j++) valShape[j]=new double[%i];")%(fs.points.size1(),fs.points.size1())
        print ("                 recurShape **children=new recurShape*[%i];")%nbChilds[iGeom]
        print ("                 for (int j=0; j<%i; j++) children[j]=NULL;")%nbChilds[iGeom]
        print ("                 for (int j=0; j<%i; j++){")%(fs.points.size1())
        print ("                    for (int k=0; k<%i; k++){")%(fs.points.size1())
        print ("                       valShape[j][k]=0;");
        print ("                       for (int l=0; l<%i; l++){")%(fs.points.size1())
        print ("                          valShape[j][k]+=childVal[i][j][l]*parent.getValShape()[l][k];")
        print ("                       }");
        print ("                    }");
        print ("                 }");
        print ("                 parent.getChildren()[i]=new recurShape(valShape,children);")
        print ("              }");
        print ("              recurChild(leftlevel,*parent.getChildren()[i]);");
        print ("          }");
        print ("       }");
        print ("   }");
        print ("   %s%i (int maxlevel){")%(geomTypes[iGeom],allN[iGeom][i])
        print ("      maxLevel=maxlevel;");
        print ("      double **valShape=new double*[%i];")%(fs.points.size1())
        print ("      for (int i=0; i<%i; i++) valShape[i]=new double[%i];")%(fs.points.size1(),fs.points.size1())
        print ("      for (int i=0; i<%i; i++) for (int j=0; j<%i; j++) if (i==j) valShape[i][j]=1; else valShape[i][j]=0;")%(fs.points.size1(),fs.points.size1())
        print ("      recurShape **children=new recurShape*[%i];")%nbChilds[iGeom]
        print ("      for (int i=0; i<%i; i++) children[i]=NULL;")%nbChilds[iGeom]
        print ("      tree=new recurShape(valShape,children);")
        print ("      recurChild(maxLevel,*tree);")
        print ("   }");
        print ("   static %s%i *_singleton;")%(geomTypes[iGeom],allN[iGeom][i])
        print ("public:")
        print ("   static %s%i *getInstance(int maxlevel){")%(geomTypes[iGeom],allN[iGeom][i])
        print ("      if (maxlevel>maxLevel && _singleton!=NULL){ //We need to extend the level")
        print ("         maxLevel=maxlevel;")
        print ("         _singleton->recurChild(maxLevel,*tree);")
        print ("      }")
        print ("      if (_singleton==NULL)")
        print ("         _singleton =  new %s%i(maxlevel);")%(geomTypes[iGeom],allN[iGeom][i])
        print ("      return _singleton;");
        print ("   }");
        print("   recurShape *getTree() {return tree;}")
        print ("};");
        print ("vtkHighOrder::shapeFunctions::%s%i *vtkHighOrder::shapeFunctions::%s%i::_singleton = NULL;")%(geomTypes[iGeom],allN[iGeom][i],geomTypes[iGeom],allN[iGeom][i])
        print ("int vtkHighOrder::shapeFunctions::%s%i::maxLevel = 0;")%(geomTypes[iGeom],allN[iGeom][i])
        print ("vtkHighOrder::recurShape *vtkHighOrder::shapeFunctions::%s%i::tree = NULL;")%(geomTypes[iGeom],allN[iGeom][i])
        print "//Initialization of static arrays"
        print outclass.getvalue(),
    print ""
    print ("vtkHighOrder::recurShape *vtkHighOrder::get%sTree(int nbpoints, int maxlevel){")%geomTypes[iGeom]
    print "   switch(nbpoints) {"
    for i in range(len(allGmshTypes[iGeom])):
        print ("   case %i:")%allN[iGeom][i]
        print ("      return vtkHighOrder::shapeFunctions::%s%i::getInstance(maxlevel)->getTree();")%(geomTypes[iGeom],allN[iGeom][i])
    print ("   default:")
    print "      vtkErrorMacro(<< \"Incompatible number of high order functions for this element\");"
    print "      return NULL;"
    print "   }"
    print "}"
