#ifndef NETWORK_EVOLVE_VISUALIZE_GUI_H
#define NETWORK_EVOLVE_VISUALIZE_GUI_H

#include <QApplication>
#include <QCoreApplication>
#include <QWidget>
#include <QMainWindow>
#include <QSettings>
#include <QLabel>
#include <QTextEdit>
#include <QPushButton>
#include <QSplitter>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QGroupBox>
#include <QComboBox>
#include <QLineEdit>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QToolBar>
#include <QFile>
#include <QProcess>
#include <QTextStream>
#include <QFileInfo>
#include <QDir>
#include <QFileDialog>
#include <QLibrary>
#include <QThread>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QGraphicsRectItem>
#include <QGraphicsLineItem>
#include <QSemaphore>
#include "CodeEditor.h"

namespace NetworkEvolutionLib
{
	class MainWindow : public QMainWindow
	{
		Q_OBJECT
		
	public:
		MainWindow();
		~MainWindow();
		QSize sizeHint() const;
	
		static int MainCallback(int, int, void**, double*, int ***);
	
	public slots:
		void go();
		void updateScene(int , int, void**, double*, int *** );
		
	private:
		
		void insertTextItem(double, int, int, int, int);
		QGraphicsScene * scene;
		QList<QGraphicsItem*> previousLine;
		
	};
	
	class Thread : public QThread
	{
			Q_OBJECT
			
		public:
			Thread(const QString&, QObject * );
			void emitSignal(int, int, void**, double*, int ***);
		signals:
			void updateScene(int , int, void**, double*, int *** );
		protected:
			QLibrary * lib;
			void run();
	};
	
	class GraphicsView : public QGraphicsView
	{
		public:
			GraphicsView(QGraphicsScene * scene = 0, QWidget * parent = 0) : QGraphicsView(scene, parent) {}
		protected:
			virtual void wheelEvent(QWheelEvent * event)
			{
				if (event->modifiers() & Qt::ControlModifier)
				{
					double factor = 1.0 + 0.2 * qAbs(event->delta()/120.0);
					if (event->delta() > 0)
						scale(factor,factor);
					else
						scale(1.0/factor,1.0/factor);
				}
				else
				{
					QGraphicsView::wheelEvent(event);
				}
			}
	};
}

#endif

