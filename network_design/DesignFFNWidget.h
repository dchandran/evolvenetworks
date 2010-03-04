#ifndef DESIGNFFNWIDGET_H
#define DESIGNFFNWIDGET_H

#include <QMainWindow>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsLineItem>
#include <QGraphicsRectItem>
#include <QGridLayout>

namespace FeedforwardNetworkDesigner
{
	class Grid : public QGraphicsScene
	{
		Q_OBJECT

	public:
		Grid(QObject * parent = 0);
		~Grid();

		void updateGridLines(const QList<double>& horizontal, const QList<double>& vertical);
		void toggleSquare(double, double, bool toggle=true);
		QRectF getRectAt(int,int);
		
	protected:
		void mousePressEvent ( QGraphicsSceneMouseEvent * mouseEvent );
		
	private:
		qreal size, gap;
		QGraphicsRectItem * boundingRect;
		
		QList<QGraphicsLineItem*> horizontalItems, verticalItems;
		QList<QGraphicsRectItem*> squares;		
	};
	
	class Table : public QWidget
	{
		Q_OBJECT

	public:
		Table();
		QSize sizeHint() const;
		
	private slots:
		void tableChanged(int,int);
		void addRow();
		void addCol();
		void removeRow();
		void removeCol();
		
	private:
		QList<double> horizontalValues, verticalValues;
	};
	
};

#endif

