#include <QDebug>
#include <QGraphicsSceneMouseEvent>
#include <QBrush>
#include <QPen>
#include "DesignFFNWidget.h"

namespace FeedforwardNetworkDesigner
{
	QList< QList<qreal> > ON_OFF_VALUES;	
	int ROWS, COLS;
	
	Grid::Grid(QObject * parent) : QGraphicsScene(parent)
	{
		ON_OFF_VALUES.clear();
		ROWS = COLS = 0;
		size = 1000.0;
		gap = 0.9;
		QRectF rect(0,0,size,size);
		
		boundingRect = new QGraphicsRectItem(rect);
		boundingRect->setPen(QPen(QBrush(QColor(100,100,255)),5));
		
		setSceneRect(rect);
		addItem(boundingRect);
	}
	
	Grid::~Grid()
	{
		if (boundingRect)
			delete boundingRect;
		for (int i=0; i < horizontalItems.size(); ++i)
			if (horizontalItems[i])
				delete horizontalItems[i];
		for (int i=0; i < verticalItems.size(); ++i)
			if (verticalItems[i])
				delete verticalItems[i];
			
		for (int i=0; i < squares.size(); ++i)
			if (squares[i])
				delete squares[i];
	}
	
	QRectF Grid::getRectAt(int m,int n)
	{
		QRectF rect;
		
		if (horizontalItems.size() > m)
			rect.
			
	}
	
	void Grid::updateGridLines(const QList<double>& horizontal, const QList<double>& vertical)
	{
		QList< QList<qreal> > temp = ON_OFF_VALUES;
		ON_OFF_VALUES.clear();
		ROWS = horizontal.size();
		COLS = vertical.size();
		
		for (int i=0; i < horizontal.size(); ++i)
		{
			ON_OFF_VALUES << QList<qreal>();
			for (int j=0; j < vertical.size(); ++j)
			{
				ON_OFF_VALUES[i] << 0.0;
			}
		}
		toggleSquare(0.0,0.0,false);
		
		
		for (int i=horizontal.size(); i < horizontalItems.size(); ++i)
			horizontalItems[i]->setVisible(false);

		for (int i=vertical.size(); i < verticalItems.size(); ++i)
			verticalItems[i]->setVisible(false);
			
		QGraphicsLineItem * line;
		
		for (int i=verticalItems.size(); i < vertical.size(); ++i)
		{
			line = new QGraphicsLineItem;
			line->setPen(boundingRect->pen());
			verticalItems << line;
			addItem(line);
		}

		for (int i=horizontalItems.size(); i < horizontal.size(); ++i)
		{
			line = new QGraphicsLineItem;
			line->setPen(boundingRect->pen());
			horizontalItems << line;
			addItem(line);
		}
		
		qreal scalex = (gap*size)/horizontal.last(), scaley = (gap*size)/vertical.last();
		qreal x,y;
		
		for (int i=0; i < horizontalItems.size() && i < horizontal.size(); ++i)
			if (horizontalItems[i])
			{	
				x = horizontal[i];
				horizontalItems[i]->setLine(0.0, x * scalex, size, x * scalex);
				horizontalItems[i]->setVisible(true);
			}
		
		for (int i=0; i < verticalItems.size() && i < vertical.size(); ++i)
			if (verticalItems[i])
			{	
				y = vertical[i];
				verticalItems[i]->setLine(y * scaley, 0.0, y * scaley, size);
				verticalItems[i]->setVisible(true);
			}
	}
	
	void Grid::toggleSquare(double x, double y, bool toggle)
	{
		if (ROWS < 1 || COLS < 1) return;
		
		QPointF p(x,y);
		int m = (int)( y / size * ROWS );
		int n = (int)( x / size * COLS );
		
		if (ROWS <= m || COLS <= n) return;
		
		if (toggle)
			qDebug() << "toggled (" << m << " , " << n << ")";

		QGraphicsRectItem * square = 0;
		
		for (int i=0; i < squares.size(); ++i)
		{
			square = squares[i];
			if (square && square->sceneBoundingRect().contains(p))
			{
				if (square->isVisible())
				{
					if (toggle)
					{
						square->setVisible(false);
						ON_OFF_VALUES[m][n] = 0.0;
					}
					else
						ON_OFF_VALUES[m][n] = 1.0;
				}
				else
				{
					if (toggle)
					{
						square->setVisible(true);
						ON_OFF_VALUES[m][n] = 1.0;
					}
					else
						ON_OFF_VALUES[m][n] = 0.0;
				}
				break;
			}
			else
				square = 0;
		}
		if (!square && toggle)
		{
			square = new QGraphicsRectItem;			
			square->setPen(Qt::NoPen);
			square->setBrush(QBrush(QColor(255,10,10)));
			square->setRect( size / COLS * n + 2.0, size / ROWS * m + 2.0, size/COLS - 8.0, size/ROWS - 8.0);
			squares << square;
			
			addItem(square);
			ON_OFF_VALUES[m][n] = 1.0;
		}
	}
	
	void Grid::mousePressEvent (QGraphicsSceneMouseEvent * mouseEvent)
	{
		QPointF pos = mouseEvent->scenePos();
		toggleSquare( pos.x(), pos.y() );
	}
	
	QSize Table::sizeHint() const
	{
		return QSize(400,400);
	}
	
	Table::Table()
	{
		Grid * grid = new Grid(this);
		QGraphicsView * view = new QGraphicsView(grid);
		view->fitInView(grid->sceneRect().adjusted(-50,-50,50,50));
		grid->updateGridLines( QList<double>() << 0.8 << 1.0,  QList<double>() << 0.5 << 1.0 );
		QGridLayout * layout = new QGridLayout;
		layout->addWidget(view);
		setLayout(layout);
	}
	
	void Table::tableChanged(int,int)
	{
	}	
	
	void Table::addRow()
	{
	}
	
	void Table::addCol()
	{
	}
	
	void Table::removeRow()
	{
	}
	
	void Table::removeCol()
	{
	}

};

