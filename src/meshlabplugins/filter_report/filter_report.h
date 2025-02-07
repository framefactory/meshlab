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

#ifndef FILTER_REPORT_H
#define FILTER_REPORT_H

#include <QJsonDocument>
#include <QJsonObject>

#include <common/interfaces.h>

class FilterReportPlugin : public MeshLabFilterInterface
{
    Q_OBJECT
    MESHLAB_PLUGIN_IID_EXPORTER(MESHLAB_FILTER_INTERFACE_IID)
    Q_INTERFACES(MeshLabFilterInterface)

public:
    bool applyFilter( const QString& filterName,MeshDocument& md,EnvWrap& env, vcg::CallBackPos * cb );

private:
    void calculateTopologicalMeasures(MeshDocument& md);
    void calculateGeometricMeasures(MeshDocument& md);
    void writeJsonToLog(const QJsonDocument& jsonDoc);
    void writeJsonToFile(const QJsonDocument& jsonDoc);

	Matrix33m computePrincipalAxisCloud(CMeshO & m);

    QJsonObject jsonReport;
    bool isWatertight;
};

#endif // FILTER_REPORT_H
