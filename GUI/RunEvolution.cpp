#include "RunEvolution.h"

using namespace NetworkEvolutionLib;
using namespace Tinkercell;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
	
	QString appDir = QCoreApplication::applicationDirPath();
	
    QFile styleFile(appDir + QString("/networkevolve.qss"));
	
	if (styleFile.open(QFile::ReadOnly | QFile::Text))
    {
        app.setStyleSheet(styleFile.readAll());
        styleFile.close();
    }

    MainWindow mainWindow;

    mainWindow.show();
    
	int output = app.exec();
	
    return output;
}

namespace NetworkEvolutionLib
{
	QSize MainWindow::sizeHint() const
	{
		return QSize(400,400);
	}

	MainWindow::MainWindow()
	{
		QGraphicsView * view = new QGraphicsView(scene = new QGraphicsScene, this);
		scene->setBackgroundBrush(QBrush(Qt::black));
		view->setDragMode(QGraphicsView::ScrollHandDrag);
		setCentralWidget(view);
		
		insertTextItem(5.0,1,5,0);
		
		insertTextItem(12.0,2,1,5);
		insertTextItem(10.0,2,2,5);
		insertTextItem(1.0,2,3,5);
		insertTextItem(42.0,2,4,5);
	}
	
	MainWindow::~MainWindow()
	{
		QList<QGraphicsItem*> items = scene->items();
		for (int i=0; i < items.size(); ++i)
			if (items[i] && !items[i]->parentItem())
				delete items[i];
	}
	
	void MainWindow::insertTextItem(double number, int line, int rank, int parent)
	{
		QGraphicsSimpleTextItem * textItem = new QGraphicsSimpleTextItem(QString::number(number));
		
		textItem->setBrush(QBrush(QColor(200,255,200)));
		textItem->setPos( 10.0 + rank * 100.0, 10.0 + line * 100.0 );
		
		QGraphicsRectItem * rectItem = new QGraphicsRectItem(	10.0 + rank * 100.0, 10.0 + line * 100.0,
																2.0 * textItem->boundingRect().width(), 
																2.0 * textItem->boundingRect().height());
		
		textItem->scale(2.0,2.0);
		
		rectItem->setBrush(QBrush(QColor(100,100,100)));
		rectItem->setPen(QPen(QColor(100,255,100)));
		
		scene->addItem(rectItem);
		scene->addItem(textItem);
		
		if (parent > 0)
		{
			QGraphicsLineItem * lineItem = new QGraphicsLineItem(
									10.0 + parent * 100.0, 30.0 + (line-1) * 100.0,
									10.0 + rank * 100.0, 10.0 + (line) * 100.0);
			lineItem->setPen(QPen(QColor(100,255,100)));
			scene->addItem(lineItem);
		}
	}
	
	void MainWindow::loadFile(const QString&)
	{
	}	
}
