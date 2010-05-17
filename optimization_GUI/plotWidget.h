#ifndef JSIMOPTIMPLOTWIDGET_H
#define JSIMOPTIMPLOTWIDGET_H

#include <QButtonGroup>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QColorDialog>
#include <QPen>
#include <QList>
#include <QColor>
#include <QDialog>
#include "DataTable.h"
#include "qwt_plot.h"
#include "qwt_color_map.h"
#include "qwt_plot_marker.h"
#include "qwt_plot_curve.h"
#include "qwt_legend.h"
#include "qwt_data.h"
#include "qwt_text.h"
#include "qwt_plot_layout.h"
#include "qwt_plot_zoomer.h"
#include "qwt_legend_item.h"
#include "qwt_scale_engine.h"

class DataPlot;

class DataColumn : public QwtData
{
public:
	DataColumn(DataTable<qreal>* data, int,int,int dt=1);
	virtual QwtData * copy() const;
	virtual size_t size() const;
	virtual double x(size_t index) const;
	virtual double y(size_t index) const;
private:
	DataTable<qreal> * dataTable;
	int column, xaxis, dt;
	
	friend class DataPlot;
};

class DataPlot : public QwtPlot
{
	Q_OBJECT
public:
	DataPlot(QWidget * parent = 0);
	void plot(const DataTable<qreal>&,int,const QString&,const QList<QPen> pens = QList<QPen>());
	virtual QSize minimumSizeHint() const;
	virtual QSize sizeHint() const;
	virtual void setLogX(bool);
	virtual void setLogY(bool);
	void replotUsingHideList();
	DataTable<qreal>& getDataTable() { return dataTable; }
	static QList<QPen> penList;
protected:
	DataTable<qreal> dataTable;
	QwtPlotZoomer * zoomer;
	QStringList hideList;
	int xcolumn, delta;
	
protected slots:
	void itemChecked(QwtPlotItem *,	bool);
	void setXAxis(int);
};


#endif
