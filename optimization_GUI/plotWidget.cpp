#include "plotWidget.h"
/**********************
	Data Column
**********************/
DataColumn::DataColumn(DataTable<qreal>* dataPtr, int xindex, int yindex, int delta)
{
	dataTable = dataPtr;
	column = yindex;
	xaxis = xindex;
	dt = delta;
	if (dt < 1 || dt > dataTable->rows()/2) dt = 1;
}

QwtData * DataColumn::copy() const
{
	return new DataColumn(dataTable,xaxis, column);
}

size_t DataColumn::size() const
{
	if (!dataTable) return 0;
	return (int)(dataTable->rows()/dt);
}

double DataColumn::x(size_t index) const
{
	if (!dataTable) return 0;
	if (xaxis < 0) return (int)index;
	return dataTable->at((int)index,xaxis);
}

double DataColumn::y(size_t index) const
{
	if (!dataTable) return 0;
	return dataTable->at((int)index*dt,column);
}

/****************************
	Data Plot
****************************/

QList<QPen> DataPlot::penList = QList<QPen>();

DataPlot::DataPlot(QWidget * parent) : QwtPlot(parent)
{
	setFrameShadow ( QFrame::Plain );
	setFrameShape ( QFrame::NoFrame );
	setFrameStyle (QFrame::NoFrame );
	xcolumn = 0;
	delta = 1;
	zoomer = new QwtPlotZoomer(xBottom,yLeft,QwtPicker::DragSelection,QwtPicker::AlwaysOff,canvas());
	zoomer->setRubberBandPen(QPen(Qt::black));
	setCanvasBackground(Qt::white);
	QwtPlotCanvas * c = canvas();
	c->setFrameShadow ( QFrame::Plain );
	c->setFrameShape ( QFrame::NoFrame );
	c->setFrameStyle (QFrame::NoFrame );
	plotLayout()->setAlignCanvasToScales(true);
	plotLayout()->setCanvasMargin(0);
	//setAxisAutoScale(xBottom);
	//setAxisAutoScale(yLeft);
	if (DataPlot::penList.isEmpty())
	{
		DataPlot::penList 	<< QPen(QColor(tr("#232CE6")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#CA420D")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#11A306")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#BF0CB0")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#D9C11F")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#0CBDBF")),2,Qt::SolidLine)
							<< QPen(QColor(tr("#232CE6")),2,Qt::DotLine)
							<< QPen(QColor(tr("#CA420D")),2,Qt::DotLine)
							<< QPen(QColor(tr("#11A306")),2,Qt::DotLine)
							<< QPen(QColor(tr("#BF0CB0")),2,Qt::DotLine)
							<< QPen(QColor(tr("#D9C11F")),2,Qt::DotLine)
							<< QPen(QColor(tr("#0CBDBF")),2,Qt::DotLine);
	}
}

QSize DataPlot::minimumSizeHint() const
{
	return QSize(100,100);
}

QSize DataPlot::sizeHint() const
{
	return QSize(160,160);
}

void DataPlot::itemChecked(QwtPlotItem * plotItem,bool on)
{
	if (plotItem)
	{
		on = !on;
		plotItem->setVisible(on);
		const QwtPlotItemList& list = itemList();
		for (int i=0; i < dataTable.cols() && i < list.size(); ++i)
			if (list.at(i) == plotItem)
			{
				if (on && hideList.contains(dataTable.colName(i)))
				{
					hideList.removeAll(dataTable.colName(i));
				}
				else
				if (!on && !hideList.contains(dataTable.colName(i)))
				{
					hideList += (dataTable.colName(i));
				}
			}
		this->replot();
	}
}

void DataPlot::setXAxis(int x)
{
	if (x >= 0 && x < dataTable.cols())
	{
		plot(dataTable,x,title().text());
	}
}

void DataPlot::plot(const DataTable<qreal>& dat, int x, const QString& title,const QList<QPen> penlist)
{
	QList<QPen> pens = penlist;

	if (pens.size() < (dat.cols()-1))
		pens = DataPlot::penList;

	delta = 1;
	xcolumn = x;
	/*if (!this->isVisible())
	{
		if (this->parentWidget() && !this->parentWidget()->isVisible())
			this->parentWidget()->show();
		else
			this->show();
	}*/
	setAutoReplot(false);
	this->dataTable = dat;

	this->clear();
	insertLegend(new QwtLegend(this), QwtPlot::RightLegend,0.2);
	legend()->setItemMode(QwtLegend::CheckableItem);
	
	connect(this,SIGNAL(legendChecked(QwtPlotItem*,bool)),this,SLOT(itemChecked(QwtPlotItem*,bool)));
	
	QList<QwtPlotCurve*> curves;
	for (int i=0, c = 0, t = 0; i < dataTable.cols(); ++i)
	{
		if (i != x && dataTable.colName(i).toLower() != tr("time"))
		{
			QwtPlotCurve * curve = new QwtPlotCurve(dataTable.colName(i));
			curve->setRenderHint(QwtPlotItem::RenderAntialiased);
			
			if (c >= pens.size())
			{
				c = 0;
			}
			
			curve->setPen(pens[c]);
			curve->setData( DataColumn(&dataTable,x,i,delta) );
			curve->attach(this);
			curve->updateLegend(legend());
			
			++c;
		}
	}
	if (dataTable.cols() > x)
		setAxisTitle(xBottom, dataTable.colName(x));
	else
		if (x < 0)
			setAxisTitle(xBottom, "Index");
		else
			setAxisTitle(xBottom, "Time");
	
	QString ylabel = axisTitle(QwtPlot::yLeft).text();
	if (ylabel.isEmpty() || ylabel.isNull())
		ylabel = tr("Values");
	setAxisTitle(yLeft, ylabel);
	setTitle(title);

	setAxisAutoScale(xBottom);
	setAxisAutoScale(yLeft);
	setAutoReplot(true);
	replot();
	if (zoomer)
	{
		zoomer->setZoomBase();
	}
	
	replotUsingHideList();
}

void DataPlot::replotUsingHideList()
{
	const QwtPlotItemList& list = itemList();
	QwtLegend * leg = legend();
	for (int i=0; i < dataTable.cols() && i < list.size(); ++i)
		if (hideList.contains(dataTable.colName(i)))
		{
			list[i]->setVisible(false);
			QWidget * w = leg->find(list[i]);
			if ( w && w->inherits( "QwtLegendItem" ) )
				((QwtLegendItem *)w)->setChecked( true );
		}
	replot();
}

void DataPlot::setLogX(bool b)
{
	if (b)
	{
		setAxisMaxMajor(QwtPlot::xBottom, 6);
		setAxisMaxMinor(QwtPlot::xBottom, 10);
		QwtLog10ScaleEngine * engine = new QwtLog10ScaleEngine;
		double d1,d2,d3;
		engine->autoScale(1,d1,d2,d3);
		setAxisScaleEngine(QwtPlot::xBottom, engine);
	}
	else
	{
		setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine);
	}
	replot();
}

void DataPlot::setLogY(bool b)
{
	if (b)
	{
		setAxisMaxMajor(QwtPlot::yLeft, 6);
		setAxisMaxMinor(QwtPlot::yLeft, 10);
		QwtLog10ScaleEngine * engine = new QwtLog10ScaleEngine;
		double d1,d2,d3;
		engine->autoScale(1,d1,d2,d3);
		setAxisScaleEngine(QwtPlot::yLeft, engine);
	}
	else
	{
		setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine);
	}
	replot();
}
	