#ifndef DESIGNFFNWIDGET_H
#define DESIGNFFNWIDGET_H

#include <QMainWindow>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsLineItem>

namespace FeedforwardNetworkDesigner
{
	class Grid : public QGraphicsScene
	{
		Q_OBJECT

	public:
		Grid();
		void updateGrid(const QList<double>& horizontal, const QList<double>& vertical);

	private:
		QList<QGraphicsLineItem*> horizontalItems, verticalItems;
		int size;
	};
	
	class Table : public QWidget
	{
		Q_OBJECT

	public:
		Table();
		
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

