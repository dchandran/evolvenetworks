#include "NetworkEvolveGUI.h"

using namespace NetworkEvolutionLib;


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
	MainWindow::MainWindow()
	{
		QSplitter * twoCols = new QSplitter;
		twoCols->setOrientation ( Qt::Horizontal );
		
		QSplitter * firstCol = new QSplitter;
		firstCol->setOrientation ( Qt::Vertical );
		
		QVBoxLayout * layout1 = new QVBoxLayout;
		QGroupBox * group1 = new QGroupBox;
		group1->setTitle(" Configure Network Types ");
		layout1->addWidget(setupNetworkOptions());
		group1->setLayout(layout1);
		
		QVBoxLayout * layout2 = new QVBoxLayout;
		QGroupBox * group2 = new QGroupBox;
		group2->setTitle(" Configure Runs ");
		layout2->addWidget(setupGAOptions());
		group2->setLayout(layout2);
		
		firstCol->addWidget(group1);
		firstCol->addWidget(group2);
		
		QVBoxLayout * layout3 = new QVBoxLayout;
		QGroupBox * group3 = new QGroupBox;
		group3->setTitle(" Fitness Function ");
		layout3->addWidget(setupEditor());
		group3->setLayout(layout3);
		
		twoCols->addWidget(firstCol);
		twoCols->setStretchFactor(0,0);
		
		twoCols->addWidget(group3);
		twoCols->setStretchFactor(1,20);
		
		setWindowTitle("Network Evolution GUI");
		
		twoCols->setStyleSheet("background-color: transparent");
		
		QToolBar * toolbar = new QToolBar(this);
		
		QPushButton * button;
		
		button = new QPushButton(toolbar);
		button->setText("RUN");
		connect(button,SIGNAL(pressed()),this,SLOT(run()));
		toolbar->addWidget(button);
		
		button = new QPushButton(toolbar);
		button->setText("Reset");
		toolbar->addWidget(button);
		
		button = new QPushButton(toolbar);
		button->setText("Quit");
		connect(button,SIGNAL(pressed()),this,SLOT(close()));
		toolbar->addWidget(button);
		
		toolbar->setMinimumHeight(40);
		
		addToolBar(Qt::BottomToolBarArea,toolbar);
		
		setCentralWidget(twoCols);
	}
	
	MainWindow::~MainWindow()
	{
	}
	
	void MainWindow::setupMassActionNetwork(QTreeWidget* treeWidget)
	{
		QTreeWidgetItem * massaction = new QTreeWidgetItem;
		massaction->setText(0,"Mass Action");
		massaction->setToolTip(0,"evolve networks with mass action kinetics");
		treeWidget->addTopLevelItem(massaction);
		
		QCheckBox * checkBox;
		treeWidget->setItemWidget(massaction,1,checkBox = new QCheckBox);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useMassAction(bool)));		
		
		QDoubleSpinBox * doubleSpinBox;
		QTreeWidgetItem * child;
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be mass action");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setMassActionProb(double)));
		doubleSpinBox->setValue(0.0);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Uni-Uni");
		child->setToolTip(0,"Proportion of the reactions that will have one reactant and product");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setUniUni(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Uni-Bi");
		child->setToolTip(0,"Proportion of the reactions that will have one reactant and two products");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setUniBi(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Bi-Uni");
		child->setToolTip(0,"Proportion of the reactions that will have two reactants and one product");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setBiUni(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Bi-Bi");
		child->setToolTip(0,"Proportion of the reactions that will have two reactants and products");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setBiBi(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent No reactants");
		child->setToolTip(0,"Proportion of the reactions that represent flux into the system (no reactant). This proportion is independent of the other proportions (uni-uni, bi-uni, etc.)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setNoReactant(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent No products");
		child->setToolTip(0,"Proportion of the reactions that represent flow out of the system (no product). This proportion is independent of the other proportions (uni-uni, bi-uni, etc.)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setNoProduct(double)));
		doubleSpinBox->setValue(0.2);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. rate constant");
		child->setToolTip(0,"Average value for a reaction rate constant");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_init_max_constant(double)));
		doubleSpinBox->setValue(1.0);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate rate constant");
		child->setToolTip(0,"Probability of mutating a rate constant during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_constants(double)));
		doubleSpinBox->setValue(0.5);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate remove reaction");
		child->setToolTip(0,"Probability of removing a reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_remove_reaction(double)));
		doubleSpinBox->setValue(0.25);
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate add reaction");
		child->setToolTip(0,"Probability of adding a new reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_add_reaction(double)));
		doubleSpinBox->setValue(0.25);
	}
	
	void MainWindow::setupEnzymeNetwork(QTreeWidget* treeWidget)
	{
		QTreeWidgetItem * enzyme = new QTreeWidgetItem; 
		enzyme->setText(0,"Enzyme Catalyzed");
		enzyme->setToolTip(0,"evolve networks with mass action and reversible Hill equation kinetics");
		treeWidget->addTopLevelItem(enzyme);
		
		QCheckBox * checkBox;
		QDoubleSpinBox * doubleSpinBox;
		QTreeWidgetItem * child;
		
		treeWidget->setItemWidget(enzyme,1,checkBox = new QCheckBox);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useEnzyme(bool)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be enzymatic networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setEnzymeProb(double)));
		doubleSpinBox->setValue(0.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Kcat range");
		child->setToolTip(0,"Initial range for the Kcat parameter");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_kcat(double)));
		doubleSpinBox->setValue(4.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Keq range (log)");
		child->setToolTip(0,"Initial range for the Keq parameter in log scale (base 2)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_log_keq(double)));
		doubleSpinBox->setValue(4.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Alpha range (log)");
		child->setToolTip(0,"Initial range for the Alpha parameter in log scale (base 2)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_alpha(double)));
		doubleSpinBox->setValue(2.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Hill range");
		child->setToolTip(0,"Initial range for the Hill coefficient");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_h);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_h(double)));
		doubleSpinBox->setValue(4.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"S-half range");
		child->setToolTip(0,"Initial range for the substrate half saturation value");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_s_half);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_s_half(double)));
		doubleSpinBox->setValue(20.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"S-half range");
		child->setToolTip(0,"Initial range for the substrate half saturation value");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_init_max_p_half(double)));
		doubleSpinBox->setValue(20.0);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate enzyme");
		child->setToolTip(0,"Probability of swapping an enzyme during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_enzyme(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Kcat");
		child->setToolTip(0,"Probability of mutating a Kcat (catalytic constant) during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_k_cat(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Keq");
		child->setToolTip(0,"Probability of mutating the equilibrium constant during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_k_eq(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Alpha");
		child->setToolTip(0,"Probability of mutating the Alpha value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_alpha(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Hill");
		child->setToolTip(0,"Probability of mutating the Hill coefficient during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_h(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate S-half");
		child->setToolTip(0,"Probability of mutating the substrate half-saturation point during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_s_half(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate P-half");
		child->setToolTip(0,"Probability of mutating the product half-saturation point during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setValue(enzyme_init_max_p_half);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_p_half(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate remove reaction");
		child->setToolTip(0,"Probability of removing a reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_remove(double)));
		doubleSpinBox->setValue(0.1);
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate add reaction");
		child->setToolTip(0,"Probability of adding a new reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_add(double)));
		doubleSpinBox->setValue(0.1);
	}
	
	void MainWindow::setupProteinNetwork(QTreeWidget* treeWidget)
	{
		QTreeWidgetItem * protein = new QTreeWidgetItem;
		protein->setText(0,"Protein State Change");
		protein->setToolTip(0,"evolve networks with proteins switching between active and inactive states");
		treeWidget->addTopLevelItem(protein);
		
		QCheckBox * checkBox;
		QDoubleSpinBox * doubleSpinBox;
		QTreeWidgetItem * child;
		
		treeWidget->setItemWidget(protein,1,checkBox = new QCheckBox);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useProtein(bool)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be protein interaction networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setProteinProb(double)));
		doubleSpinBox->setValue(0.0);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. Km");
		child->setToolTip(0,"The average Michaelis-Menten constant for initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_init_ka(double)));
		doubleSpinBox->setValue(10.0);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. Vmax");
		child->setToolTip(0,"The average Vmax constant for initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_init_vmax(double)));
		doubleSpinBox->setValue(10.0);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. total");
		child->setToolTip(0,"The average total concentration for proteins in the initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_init_total(double)));
		doubleSpinBox->setValue(10.0);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate regulation");
		child->setToolTip(0,"Probability of rewiring the network during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_rewire(double)));
		doubleSpinBox->setValue(0.2);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate parameter");
		child->setToolTip(0,"Probability of changing a parameter during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_parameter(double)));
		doubleSpinBox->setValue(0.4);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate total");
		child->setToolTip(0,"Probability of changing a conservation law (total concentration) during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_total(double)));
		doubleSpinBox->setValue(0.2);
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Add/remove proteins");
		child->setToolTip(0,"Probability of adding or removing proteins in the network during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_addremove(double)));
		doubleSpinBox->setValue(0.2);
	}
	
	void MainWindow::setupGeneticNetwork(QTreeWidget* treeWidget)
	{
		QTreeWidgetItem * grn = new QTreeWidgetItem;
		grn->setText(0,"Gene Regulation");
		grn->setToolTip(0,"evolve gene gegulatory networks that use fractional saturation model");
		treeWidget->addTopLevelItem(grn);
		
		QCheckBox * checkBox;
		QDoubleSpinBox * doubleSpinBox;
		QSpinBox * intSpinBox;
		QTreeWidgetItem * child;
		
		treeWidget->setItemWidget(grn,1,checkBox = new QCheckBox);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useGRN(bool)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be genetic networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setGRNProb(double)));
		doubleSpinBox->setValue(1.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Resouce inflow");
		child->setToolTip(0,"Use a non-zero value here to model resource consumption. This value determines the constant inflow of resources from the environment.");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_inflow(double)));
		set_grn_init_inflow(0.0);
		doubleSpinBox->setValue(0.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Resouce consumption");
		child->setToolTip(0,"Use a non-zero value here to model resource consumption. This value determines the amount of resources consumed due to protein production.");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_cost_per_protein(double)));
		set_grn_init_cost_per_protein(0.0);
		doubleSpinBox->setValue(0.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Max. complex size");
		child->setToolTip(0,"The maximum number of transcription factors that can form a complex");
		treeWidget->setItemWidget(child,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,10);
		intSpinBox->setSingleStep(1);		
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_grn_init_max_complex_size(int)));
		intSpinBox->setValue(4);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. Ka");
		child->setToolTip(0,"Average association constant for transcription factors");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_Ka(double)));
		doubleSpinBox->setValue(10.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. Vmax");
		child->setToolTip(0,"Average Vmax for protein production rate");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_Vmax(double)));
		doubleSpinBox->setValue(10.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Avg. degradation");
		child->setToolTip(0,"Average degradation/dilution rate for proteins");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_degradation(double)));
		doubleSpinBox->setValue(2.0);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Ka");
		child->setToolTip(0,"Probability of changing a Ka value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_Ka(double)));
		doubleSpinBox->setValue(0.2);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Vmax");
		child->setToolTip(0,"Probability of changing a Vmax value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_Vmax(double)));
		doubleSpinBox->setValue(0.2);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Complex");
		child->setToolTip(0,"Probability of changing the proteins inside a complex during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_complex(double)));
		doubleSpinBox->setValue(0.2);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Add genes");
		child->setToolTip(0,"Probability of adding a new gene during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_add_gene(double)));
		doubleSpinBox->setValue(0.2);
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Remove genes");
		child->setToolTip(0,"Probability of removing a gene during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_remove_gene(double)));
		doubleSpinBox->setValue(0.2);
	}
	
	
	QWidget * MainWindow::setupNetworkOptions()
	{
		QTreeWidget * treeWidget = new QTreeWidget();
		treeWidget->setColumnCount(2);
		treeWidget->setColumnWidth(0,200);
		treeWidget->setHeaderLabels( QStringList() << "network/property" << "value" );
		
		setupMassActionNetwork(treeWidget);
		setupEnzymeNetwork(treeWidget);
		setupProteinNetwork(treeWidget);
		setupGeneticNetwork(treeWidget);
		
		
		QCheckBox * checkBox;
		QDoubleSpinBox * doubleSpinBox;
		QSpinBox * intSpinBox;
		QTreeWidgetItem * option;
		
		option = new QTreeWidgetItem;
		option->setText(0,"Avg. Molecular Species");
		option->setToolTip(0,"Average number of molecular species is the initial set of networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);		
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_species(int)));
		intSpinBox->setValue(8);
		
		option = new QTreeWidgetItem;
		option->setText(0,"Avg. Num. Reactions");
		option->setToolTip(0,"Average number of reactions is the initial set of networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);		
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_reactions(int)));
		intSpinBox->setValue(8);
		
		option = new QTreeWidgetItem;
		option->setText(0,"Crossover Rate");
		option->setToolTip(0,"Probability of a crossover event when generating the next generation of networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_crossover_rate(double)));
		doubleSpinBox->setValue(1.0);
		
		option = new QTreeWidgetItem;
		option->setText(0,"Avg. initial values");
		option->setToolTip(0,"Average concentration values of molecules in the initial networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_init_iv(double)));
		doubleSpinBox->setValue(1.0);
		
		option = new QTreeWidgetItem;
		option->setText(0,"Mutate initial values");
		option->setToolTip(0,"Probability of mutating the concentration values during a mutation event");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_mutate_iv(double)));
		doubleSpinBox->setValue(0.05);
		
		
		option = new QTreeWidgetItem;
		option->setText(0,"Lineage tracking");
		option->setToolTip(0,"Track the original parent(s) for each network, and record this table in the log file");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,checkBox = new QCheckBox);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(set_lineageTracking(bool)));
		checkBox->setChecked(true);
		
		return treeWidget;
	}
	
	QWidget * MainWindow::setupGAOptions()
	{
		QTreeWidget * treeWidget = new QTreeWidget();
		treeWidget->setHeaderLabels( QStringList() << "parameter" << "value" );
		
		treeWidget->setColumnCount(2);
		treeWidget->setColumnWidth(0,160);
		QDoubleSpinBox * doubleSpinBox;
		QSpinBox * intSpinBox;
		QLineEdit * lineEdit;
		
		QTreeWidgetItem * code = new QTreeWidgetItem;
		code->setText(0,"Code");
		code->setToolTip(0,"The code that is automatically generated");
		treeWidget->addTopLevelItem(code);
		treeWidget->setItemWidget(code,1,lineEdit = new QLineEdit(codeFile = tr("run_net_ga.c")));
		connect(lineEdit,SIGNAL(textEdited (const QString &)),this,SLOT(setCodeFile(const QString&)));
		
		QTreeWidgetItem * output = new QTreeWidgetItem;
		output->setText(0,"Log file");
		output->setToolTip(0,"This file that will contain the results from all the runs, which includes the network, fitness values, fitness distrubtion, and lineage tracking");
		treeWidget->addTopLevelItem(output);
		treeWidget->setItemWidget(output,1,lineEdit = new QLineEdit(logFile = tr("evolution.log")));
		connect(lineEdit,SIGNAL(textEdited (const QString &)),this,SLOT(setLogFile(const QString&)));
		
		QTreeWidgetItem * runs = new QTreeWidgetItem;
		runs->setText(0,"Runs");
		runs->setToolTip(0,"The number of times to repeat the evolution experiment. Each result will be different, unless the same seed is used");
		treeWidget->addTopLevelItem(runs);
		treeWidget->setItemWidget(runs,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setRuns(int)));
		intSpinBox->setValue(10);
		
		QTreeWidgetItem * popSz = new QTreeWidgetItem;
		popSz->setText(0,"Population Size");
		popSz->setToolTip(0,"The initial population size for each run");
		treeWidget->addTopLevelItem(popSz);
		treeWidget->setItemWidget(popSz,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setPopSz(int)));
		intSpinBox->setValue(200);
		
		QTreeWidgetItem * generation = new QTreeWidgetItem;
		generation->setText(0,"Generations");
		generation->setToolTip(0,"The number of iterations to run the genetic algorithm during each evolution experiment");
		treeWidget->addTopLevelItem(generation);
		treeWidget->setItemWidget(generation,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setGenerations(int)));
		intSpinBox->setValue(50);
		
		QTreeWidgetItem * seed = new QTreeWidgetItem;
		seed->setText(0,"Seed");
		seed->setToolTip(0,"Set the random number generator seed. This is used for repeating the same experiment.");
		treeWidget->addTopLevelItem(seed);
		treeWidget->setItemWidget(seed,1,lineEdit = new QLineEdit(""));
		connect(lineEdit,SIGNAL(textEdited (const QString &)),this,SLOT(setSeed(const QString&)));
		
		return treeWidget;
	}
	
	QWidget * MainWindow::setupEditor()
	{
		QComboBox * comboBox = new QComboBox;
		comboBox->addItems( QStringList() << "oscillation" << "noise" );
		
		codeEditor = new Tinkercell::CodeEditor;
		codeEditor->setStyleSheet("background-color: #F9F7E4");
		
		CSyntaxHighlighter * syntaxHighlighter = new CSyntaxHighlighter(codeEditor->document());
		
		QVBoxLayout * layout = new QVBoxLayout;
		layout->addWidget(comboBox);
		layout->addWidget(codeEditor);
		
		QWidget * widget = new QWidget;
		widget->setLayout(layout);
		
		return widget;
	}
	
	QSize MainWindow::sizeHint() const
	{
		return QSize(800,600);
	}
	
	QString MainWindow::init()
	{	
		QString s =	
				tr("void init()\n")
				+ tr("{\n ")
				+ tr("   setDistributionOfMassActionNetwork(")
				+ QString::number(uni_uni) + tr(",")
				+ QString::number(uni_bi) + tr(",")
				+ QString::number(bi_uni) + tr(",")
				+ QString::number(bi_bi) + tr(",")
				+ QString::number(no_reactant) + tr(",")
				+ QString::number(no_product)
				+ tr(");\n")
				+ tr("    setRateConstantForMassActionNetwork(")
				+ QString::number(ma_init_max_constant)
				+ tr(");\n")
				+ tr("    setMutationRatesForMassActionNetwork(")
				+ QString::number(ma_mutate_constants) + tr(",")
				+ QString::number(ma_mutate_remove_reaction) + tr(",")
				+ QString::number(ma_mutate_add_reaction)
				+ tr(");\n")
				+ tr("    setRateConstantsForProteinInteractionNetwork(")
				+ QString::number(prot_init_ka) + tr(",")
				+ QString::number(prot_init_vmax) + tr(",")
				+ QString::number(prot_init_total)
				+ tr(");\n")
				+ tr("    setMutationRatesForProteinInteractionNetwork(")
				+ QString::number(prot_mutate_rewire) + tr(",")
				+ QString::number(prot_mutate_parameter) + tr(",")
				+ QString::number(prot_mutate_total) + tr(",")
				+ QString::number(prot_mutate_addremove)
				+ tr(");\n")
				+ tr("    setRateConstantsForEnzymeNetwork(")
				+ QString::number(enzyme_init_max_kcat) + tr(",")
				+ QString::number(enzyme_init_max_log_keq) + tr(",")
				+ QString::number(enzyme_init_max_alpha) + tr(",")
				+ QString::number(enzyme_init_max_h) + tr(",")
				+ QString::number(enzyme_init_max_s_half) + tr(",")
				+ QString::number(enzyme_init_max_p_half)
				+ tr(");\n")
				+ tr("    setMutationRatesForEnzymeNetwork(")
				+ QString::number(enzyme_mutate_enzyme) + tr(",")
				+ QString::number(enzyme_mutate_k_cat) + tr(",")
				+ QString::number(enzyme_mutate_k_eq) + tr(",")
				+ QString::number(enzyme_mutate_alpha) + tr(",")
				+ QString::number(enzyme_mutate_h) + tr(",")
				+ QString::number(enzyme_mutate_s_half) + tr(",")
				+ QString::number(enzyme_mutate_p_half) + tr(",")
				+ QString::number(enzyme_mutate_remove) + tr(",")
				+ QString::number(enzyme_mutate_add)
				+ tr(");\n")
				+ tr("    setResourceRestriction(")
				+ QString::number(grn_init_inflow) + tr(",")
				+ QString::number(grn_init_cost_per_protein)
				+ tr(");\n")
				+ tr("    setRateConstantsForGeneRegulationNetwork(")
				+ QString::number(grn_init_max_complex_size) + tr(",")
				+ QString::number(grn_init_Ka) + tr(",")
				+ QString::number(grn_init_Vmax) + tr(",")
				+ QString::number(grn_init_degradation)
				+ tr(");\n")
				+ tr("    setMutationRatesForGeneRegulationNetwork(")
				+ QString::number(grn_mutate_Ka) + tr(",")
				+ QString::number(grn_mutate_Vmax) + tr(",")
				+ QString::number(grn_mutate_complex) + tr(",")
				+ QString::number(grn_mutate_add_gene) + tr(",")
				+ QString::number(grn_mutate_remove_gene)
				+ tr(");\n")
				+ tr("    setInitialNetworkSize(")
				+ QString::number(species) + tr(",")
				+ QString::number(reactions) 
				+ tr(");\n")
				+ tr("    setCrossoverRate(")
				+ QString::number(crossover_rate) 
				+ tr(");\n")
				+ tr("    setAverageInitialValue(")
				+ QString::number(init_iv) 
				+ tr(");\n")
				+ tr("    setMutationRateOfInitialValues(")
				+ QString::number(mutate_iv) 
				+ tr(");\n");
		
		if (lineageTracking)
			s += tr("    lineageTrackingON();\n");
		else
			s += tr("    lineageTrackingOFF();\n");
		
		s += tr("}\n");
		
		return s;
	}
	
	QString MainWindow::callbackFunction()
	{
		QString s = tr("int callback(int iter, GAPopulation P, int popSz)\n")
			+ tr("{\n")
			+ tr("    return 0;\n}\n");
		
		return s;
	}
	
	QString MainWindow::mainFunction()
	{
		QString s = tr("int main()\n")
			+ tr("{\n")
			+ tr("    int runs = ") + QString::number(runs) + tr(";\n")
			+ tr("    int genations = ") + QString::number(generations) + tr(";\n")
			+ tr("    int popSz = ") + QString::number(popSz) + tr(";\n")
			+ tr("    GAPopulation P;\n")
			+ tr("    int i,j;\n\n")
			+ tr("    init();\n")
			+ tr("    for (i=0; i < runs; ++i)\n")
			+ tr("       P = evolveNetworks(popSz*5,popSz,generations,&callback);\n")
			+ tr("    return 0;\n}\n");
		
		return s;
	}
	
	void MainWindow::run()
	{
		QString appDir = QCoreApplication::applicationDirPath();
		QFile qfile(appDir + tr("/") + codeFile);
		if (!qfile.open(QIODevice::WriteOnly | QIODevice::Text))
			return;

		QTextStream out(&qfile);
		out << 
			init() + tr("\n") + 
			codeEditor->toPlainText() + tr("\n") + 
			callbackFunction() + tr("\n") + 
			mainFunction();

		qfile.close();
	}
	
	void MainWindow::reset()
	{
	}
	
	void MainWindow::quit()
	{
	}
	
}
