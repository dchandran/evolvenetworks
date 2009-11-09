#include "ModuleSelector.h"

namespace NetworkEvolutionLib
{
	ModuleSelector::ModuleSelector(QWidget * parent): QWidget(parent)
	{
		QTabWidget * tabWidget = new QTabWidget;
		QHBoxLayout * layout = new QHBoxLayout;
		layout->addWidget(tabWidget);
		setLayout(layout);
	}
}
 