/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <QJsonDocument>
#include <QJsonObject>
#include <QFile>
#include <QTextStream>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/inertia.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/update/selection.h>
#include<vcg/complex/append.h>
#include<vcg/simplex/face/pos.h>
#include<vcg/complex/algorithms/bitquad_support.h>
#include<vcg/complex/algorithms/mesh_to_matrix.h>
#include<vcg/complex/algorithms/bitquad_optimization.h>

#include "filter_report.h"

using namespace std;
using namespace vcg;

// Core Function doing the actual mesh processing.
bool FilterReportPlugin::applyFilter( const QString& filterName,MeshDocument& md,EnvWrap& env, vcg::CallBackPos * /*cb*/ )
{
    if (filterName == "Generate JSON Report")
    {
        calculateTopologicalMeasures(md);
        calculateGeometricMeasures(md);

        QJsonDocument jsonDoc;
        jsonDoc.setObject(jsonReport);

        writeJsonToFile(jsonDoc);
        writeJsonToLog(jsonDoc);

        return true;
    }

    return false;
}

void FilterReportPlugin::calculateTopologicalMeasures(MeshDocument& md)
{
    CMeshO &m = md.mm()->cm;
    tri::Allocator<CMeshO>::CompactFaceVector(m);
    tri::Allocator<CMeshO>::CompactVertexVector(m);
    md.mm()->updateDataMask(MeshModel::MM_FACEFACETOPO);
    md.mm()->updateDataMask(MeshModel::MM_VERTFACETOPO);

    int edgeNonManifFFNum = tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m, true);
    int faceEdgeManif = tri::UpdateSelection<CMeshO>::FaceCount(m);
    tri::UpdateSelection<CMeshO>::VertexClear(m);
    tri::UpdateSelection<CMeshO>::FaceClear(m);

    int vertManifNum = tri::Clean<CMeshO>::CountNonManifoldVertexFF(m, true);
    tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(m);
    int faceVertManif = tri::UpdateSelection<CMeshO>::FaceCount(m);

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
    tri::Clean<CMeshO>::CountEdgeNum(m, edgeNum, edgeBorderNum, edgeNonManifNum);

    assert(edgeNonManifFFNum == edgeNonManifNum);

    bool isWatertight = edgeBorderNum == 0 && edgeNonManifNum == 0;
    int unrefVertNum = tri::Clean<CMeshO>::CountUnreferencedVertex(m);
    int connectedComponentsNum = tri::Clean<CMeshO>::CountConnectedComponents(m);

    int holeNum = 0;
    int genus = 0;
    bool isTwoManifold = edgeNonManifFFNum == 0 && vertManifNum == 0;

    // For Manifold meshes compute some other stuff
    if (isTwoManifold) {
        holeNum = tri::Clean<CMeshO>::CountHoles(m);
        genus = tri::Clean<CMeshO>::MeshGenus(m.vn - unrefVertNum, edgeNum, m.fn, holeNum, connectedComponentsNum);
    }

    // Write results to JSON data structure

    QJsonObject jsonStatistics;
    jsonStatistics.insert("vertices", m.vn);
    jsonStatistics.insert("edges", edgeNum);
    jsonStatistics.insert("faces", m.fn);

    QJsonObject jsonHealth;
    jsonHealth.insert("unreferencedVertices", unrefVertNum);

    QJsonObject jsonTopology;
    jsonTopology.insert("boundaryEdges", edgeBorderNum);
    jsonTopology.insert("connectedComponentCount", connectedComponentsNum);
    jsonTopology.insert("isTwoManifold", isTwoManifold);
    jsonTopology.insert("isWatertight", isWatertight);

    if (isTwoManifold) {
        jsonTopology.insert("holeCount", holeNum);
        jsonTopology.insert("genus", genus);
    }
    else {
        jsonTopology.insert("nonTwoManifoldEdgeCount", edgeNonManifFFNum);
        jsonTopology.insert("nonTwoManifoldVertexCount", vertManifNum);
        jsonTopology.insert("edgeIncidentFaces", faceEdgeManif);
        jsonTopology.insert("vertexIncidentFaces", faceVertManif);
    }

    jsonReport.insert("statistics", jsonStatistics);
    jsonReport.insert("health", jsonHealth);
    jsonReport.insert("topology", jsonTopology);
}

void FilterReportPlugin::calculateGeometricMeasures(MeshDocument& md)
{
    CMeshO &m = md.mm()->cm;
    QJsonObject jsonGeometry;

    // the mesh has to be correctly transformed
    if (m.Tr != Matrix44m::Identity()) {
        tri::UpdatePosition<CMeshO>::Matrix(m, m.Tr, true);
    }

    QJsonObject jsonBBSize;
    jsonBBSize.insert("x", m.bbox.DimX());
    jsonBBSize.insert("y", m.bbox.DimY());
    jsonBBSize.insert("z", m.bbox.DimZ());

    QJsonObject jsonBBMin;
    jsonBBMin.insert("x", m.bbox.min[0]);
    jsonBBMin.insert("y", m.bbox.min[1]);
    jsonBBMin.insert("z", m.bbox.min[2]);

    QJsonObject jsonBBMax;
    jsonBBMax.insert("x", m.bbox.max[0]);
    jsonBBMax.insert("y", m.bbox.max[1]);
    jsonBBMax.insert("z", m.bbox.max[2]);

    QJsonObject jsonBox;
    jsonBox.insert("size", jsonBBSize);
    jsonBox.insert("min", jsonBBMin);
    jsonBox.insert("max", jsonBBMax);

    jsonGeometry.insert("boundingBox", jsonBox);

    float area = tri::Stat<CMeshO>::ComputeMeshArea(m);

    Point3m shellBarycenter = tri::Stat<CMeshO>::ComputeShellBarycenter(m);
    QJsonObject jsonShellBC;
    jsonShellBC.insert("x", shellBarycenter[0]);
    jsonShellBC.insert("y", shellBarycenter[1]);
    jsonShellBC.insert("z", shellBarycenter[2]);

    Point3m cloudBarycenter = tri::Stat<CMeshO>::ComputeCloudBarycenter(m, false);
    QJsonObject jsonCloudBC;
    jsonCloudBC.insert("x", cloudBarycenter[0]);
    jsonCloudBC.insert("y", cloudBarycenter[0]);
    jsonCloudBC.insert("z", cloudBarycenter[0]);

    jsonGeometry.insert("shellBarycenter", jsonShellBC);
    jsonGeometry.insert("cloudBarycenter", jsonCloudBC);

    // if mesh contains a point cloud, terminate here
    if ((m.fn == 0) && (m.vn != 0)) {
        return;
    }

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
    tri::Clean<CMeshO>::CountEdgeNum(m, edgeNum, edgeBorderNum, edgeNonManifNum);
    bool isWatertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);

    Matrix33m matPCA;
    Point3m pcav;

    if (isWatertight) {
        tri::Inertia<CMeshO> inertia(m);

        // volume
        float volume = inertia.Mass();
        jsonGeometry.insert("volume", volume);

        // center of mass
        Point3m centerOfMass = inertia.CenterOfMass();
        QJsonObject jsonCenter;
        jsonCenter.insert("x", centerOfMass[0]);
        jsonCenter.insert("y", centerOfMass[1]);
        jsonCenter.insert("z", centerOfMass[2]);
        jsonGeometry.insert("centerOfMass", jsonCenter);

        // principal axis
        inertia.InertiaTensorEigen(matPCA, pcav);
    }
    else
    {
        matPCA = computePrincipalAxisCloud(m);
    }

    // write PCA axes
    QJsonArray jsonAxesArray;
    for (int i = 0; i < 3; ++i) {
        QJsonObject jsonAxis;
        jsonAxis.insert("x", matPCA[i][0]);
        jsonAxis.insert("y", matPCA[i][1]);
        jsonAxis.insert("z", matPCA[i][2]);
        jsonAxesArray.push_back(jsonAxis);
    }
    jsonGeometry.insert("principalAxes", jsonAxesArray);

    // transform mesh back to its original position
    if (m.Tr != Matrix44m::Identity()) {
        tri::UpdatePosition<CMeshO>::Matrix(m, Inverse(m.Tr), true);
    }
}

void FilterReportPlugin::writeJsonToLog(const QJsonDocument& jsonDoc)
{
    Log("@@JSONREPORT@@");
    Log(jsonDoc.toJson(QJsonDocument::Compact).toStdString().c_str());
}

void FilterReportPlugin::writeJsonToFile(const QJsonDocument& jsonDoc)
{
    QFile reportFile("mesh-topology-info.json");
    if (!reportFile.open(QIODevice::WriteOnly)) {
        return;
    }

    QTextStream outStream(&reportFile);
    outStream << jsonDoc.toJson(QJsonDocument::Indented);
    outStream.flush();
    reportFile.close();
}

// function to calculate principal axis for pointclouds or non-watertight meshes
Matrix33m FilterReportPlugin::computePrincipalAxisCloud(CMeshO & m)
{
	Matrix33m cov;
	Point3m bp(0, 0, 0);
	vector<Point3m> PtVec;
	for (CMeshO::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
	if (!(*vi).IsD()) 
	{
		PtVec.push_back((*vi).cP());
		bp += (*vi).cP();
	}

	bp /= m.vn;

	cov.Covariance(PtVec, bp);

	Matrix33m eigenvecMatrix;
	Eigen::Matrix3d em;
	cov.ToEigenMatrix(em);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(em);
	Eigen::Matrix3d c_vec = eig.eigenvectors();
	eigenvecMatrix.FromEigenMatrix(c_vec);

	return eigenvecMatrix;
}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterReportPlugin)
