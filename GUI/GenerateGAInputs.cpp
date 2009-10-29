#include "GenerateGAInputs.h"

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
	
	mainWindow.proc.terminate();
	
    return output;
}

namespace NetworkEvolutionLib
{
	MainWindow::MainWindow()
	{
		reset();
		
		QSplitter * twoCols = new QSplitter;
		twoCols->setOrientation ( Qt::Horizontal );
		
		QSplitter * firstCol = new QSplitter;
		firstCol->setOrientation ( Qt::Vertical );
		
		QVBoxLayout * layout1 = new QVBoxLayout;
		QGroupBox * group1 = new QGroupBox;
		group1->setTitle(" Configure Network Types ");
		layout1->addWidget(setupNetworkOptions());
		group1->setLayout(layout1);
		
		firstCol->addWidget(group1);
		
		QVBoxLayout * layout2 = new QVBoxLayout;
		QGroupBox * group2 = new QGroupBox;
		group2->setTitle(" Configure Runs ");
		layout2->addWidget(setupGAOptions());
		group2->setLayout(layout2);
		
		//QVBoxLayout * layoutC = new QVBoxLayout;
		//QGroupBox * groupC = new QGroupBox;
		//groupC->setTitle(" Compile command ");
		
		//QString appDir = QCoreApplication::applicationDirPath();
		//QLineEdit * compile = new QLineEdit(compileCommand);
		//layoutC->addWidget(compile);
		
		//groupC->setLayout(layoutC);
		
		QSplitter * splitter2 = new QSplitter;
		splitter2->setOrientation ( Qt::Vertical );
		splitter2->addWidget(group2);
		//splitter2->addWidget(groupC);
		
		firstCol->addWidget(splitter2);
		QPushButton * button;
		
		QVBoxLayout * layout3 = new QVBoxLayout;
		QGroupBox * group3 = new QGroupBox;
		group3->setTitle(" Fitness Function ");
		layout3->addWidget(setupEditor());
		
		QHBoxLayout * saveclearLayout = new QHBoxLayout;
		button = new QPushButton;
		button->setText("Save");
		button->setIcon(QIcon(":/save.png"));
		button->setMaximumWidth(100);
		connect(button,SIGNAL(pressed()),this,SLOT(save()));
		saveclearLayout->addWidget(button, Qt::AlignLeft);
		
		button = new QPushButton;
		button->setText("New");
		button->setIcon(QIcon(":/new.png"));
		button->setMaximumWidth(100);
		connect(button,SIGNAL(pressed()),this,SLOT(clear()));
		saveclearLayout->addWidget(button, Qt::AlignLeft);
		
		layout3->addLayout(saveclearLayout);
		group3->setLayout(layout3);
		
		twoCols->addWidget(firstCol);
		twoCols->setStretchFactor(0,1);
		
		twoCols->addWidget(group3);
		twoCols->setStretchFactor(1,2);
		
		setWindowTitle("Network Evolution GUI");
		
		twoCols->setStyleSheet("background-color: transparent");
		
		QToolBar * toolbar = new QToolBar(this);
		
		button = new QPushButton(toolbar);
		button->setText("RUN");
		button->setIcon(QIcon(":/play.png"));
		connect(button,SIGNAL(pressed()),this,SLOT(run()));
		toolbar->addWidget(button);
		
		/*button = new QPushButton(toolbar);
		button->setText("STOP");
		button->setIcon(QIcon(":/stop.png"));
		connect(button,SIGNAL(pressed()),&proc,SLOT(terminate()));
		toolbar->addWidget(button);
		*/
		
		button = new QPushButton(toolbar);
		button->setText("Quit");
		button->setIcon(QIcon(":/exit.png"));
		connect(button,SIGNAL(pressed()),this,SLOT(reset()));
		connect(button,SIGNAL(pressed()),this,SLOT(close()));
		toolbar->addWidget(button);
		
		button = new QPushButton(toolbar);
		button->setText("Save and Quit ");
		button->setIcon(QIcon(":/exit.png"));
		connect(button,SIGNAL(pressed()),this,SLOT(close()));
		toolbar->addWidget(button);
		
		toolbar->setMinimumHeight(40);
		
		addToolBar(Qt::BottomToolBarArea,toolbar);
		
		setCentralWidget(twoCols);
		
		clear();
	}
	
	MainWindow::~MainWindow()
	{
		QSettings settings("UWashington","NetworkEvolutionLib");
		settings.beginGroup("parameters");
		
		settings.setValue("uni_uni",uni_uni);
		settings.setValue("uni_bi",uni_bi);
		settings.setValue("bi_uni",bi_uni);
		settings.setValue("bi_bi",bi_bi);
		settings.setValue("no_reactant",no_reactant);
		settings.setValue("no_product",no_product);
		settings.setValue("ma_min_constant",ma_min_constant);
		settings.setValue("ma_max_constant",ma_max_constant);
		settings.setValue("ma_mutate_constants",ma_mutate_constants);
		settings.setValue("ma_mutate_remove_reaction",ma_mutate_remove_reaction);
		settings.setValue("ma_mutate_add_reaction",ma_mutate_add_reaction);
		
		settings.setValue("prot_min_ka",prot_min_ka);
		settings.setValue("prot_min_vmax",prot_min_vmax);
		settings.setValue("prot_min_total",prot_min_total);
		settings.setValue("prot_max_ka",prot_max_ka);
		settings.setValue("prot_max_vmax",prot_max_vmax);
		settings.setValue("prot_max_total",prot_max_total);
		settings.setValue("prot_mutate_rewire",prot_mutate_rewire);
		settings.setValue("prot_mutate_parameter",prot_mutate_parameter);
		settings.setValue("prot_mutate_total",prot_mutate_total);
		settings.setValue("prot_mutate_addremove",prot_mutate_addremove);
		
		settings.setValue("enzyme_min_kcat",enzyme_min_kcat);
		settings.setValue("enzyme_min_log_keq",enzyme_min_log_keq);
		settings.setValue("enzyme_min_alpha",enzyme_min_alpha);
		settings.setValue("enzyme_min_h",enzyme_min_h);
		settings.setValue("enzyme_min_s_half",enzyme_min_s_half);
		settings.setValue("enzyme_min_p_half",enzyme_min_p_half);
		settings.setValue("enzyme_max_kcat",enzyme_max_kcat);
		settings.setValue("enzyme_max_log_keq",enzyme_max_log_keq);
		settings.setValue("enzyme_max_alpha",enzyme_max_alpha);
		settings.setValue("enzyme_max_h",enzyme_max_h);
		settings.setValue("enzyme_max_s_half",enzyme_max_s_half);
		settings.setValue("enzyme_max_p_half",enzyme_max_p_half);
		
		settings.setValue("enzyme_mutate_enzyme",enzyme_mutate_enzyme);
		settings.setValue("enzyme_mutate_k_cat",enzyme_mutate_k_cat);
		settings.setValue("enzyme_mutate_k_eq",enzyme_mutate_k_eq);
		settings.setValue("enzyme_mutate_h",enzyme_mutate_h);
		settings.setValue("enzyme_mutate_s_half",enzyme_mutate_s_half);
		settings.setValue("enzyme_mutate_p_half",enzyme_mutate_p_half);
		settings.setValue("enzyme_mutate_remove",enzyme_mutate_remove);
		settings.setValue("enzyme_mutate_add",enzyme_mutate_add);
		
		settings.setValue("grn_init_inflow",grn_init_inflow);
		settings.setValue("grn_init_cost_per_protein",grn_init_cost_per_protein);
		settings.setValue("grn_min_complex_size",grn_min_complex_size);
		settings.setValue("grn_min_Ka",grn_min_Ka);
		settings.setValue("grn_min_Vmax",grn_min_Vmax);
		settings.setValue("grn_min_degradation",grn_min_degradation);
		settings.setValue("grn_max_complex_size",grn_max_complex_size);
		settings.setValue("grn_max_Ka",grn_max_Ka);
		settings.setValue("grn_max_Vmax",grn_max_Vmax);
		settings.setValue("grn_max_degradation",grn_max_degradation);
		
		settings.setValue("grn_mutate_Ka",grn_mutate_Ka);
		settings.setValue("grn_mutate_Vmax",grn_mutate_Vmax);
		settings.setValue("grn_mutate_complex",grn_mutate_complex);
		settings.setValue("grn_mutate_add_gene",grn_mutate_add_gene);
		settings.setValue("grn_mutate_remove_gene",grn_mutate_remove_gene);
		
		settings.setValue("species_min",species_min);
		settings.setValue("species_max",species_max);
		settings.setValue("reactions_min",reactions_min);
		settings.setValue("reactions_max",reactions_max);
		settings.setValue("crossover_rate",crossover_rate);
		settings.setValue("init_iv",init_iv);
		settings.setValue("mutate_iv",mutate_iv);
		settings.setValue("lineageTracking",lineageTracking);
		
		settings.setValue("mass_action_prob",mass_action_prob);
		settings.setValue("enzyme_prob",enzyme_prob);
		settings.setValue("protein_net_prob",protein_net_prob);
		settings.setValue("grn_prob",grn_prob);
		
		settings.setValue("codeFile",codeFile);
		settings.setValue("logFile",logFile);
		settings.setValue("runs",runs);
		settings.setValue("generations",generations);
		settings.setValue("popSz",popSz);
		settings.setValue("initPopSz",initPopSz);
		settings.setValue("compileCommand",compileCommand);
		
		settings.setValue("bestNetworkFitness1",bestNetworkFitness1);
		settings.setValue("bestNetworkScript1",bestNetworkScript1);
		settings.setValue("bestNetworkSize1",bestNetworkSize1);
		settings.setValue("bestNetworkLineage1",bestNetworkLineage1);
		settings.setValue("allFitness1",allFitness1);
		settings.setValue("allNetworkLineage1",allNetworkLineage1);
		
		settings.setValue("bestNetworkFitness2",bestNetworkFitness2);
		settings.setValue("bestNetworkScript2",bestNetworkScript2);
		settings.setValue("bestNetworkSize2",bestNetworkSize2);
		settings.setValue("bestNetworkLineage2",bestNetworkLineage2);
		settings.setValue("allFitness2",allFitness2);
		settings.setValue("allNetworkLineage2",allNetworkLineage2);
		settings.setValue("seeds2",seeds2);
		
		settings.setValue("max_fitness",max_fitness);
		
		settings.endGroup();
	}
	
	void MainWindow::setupMassActionNetwork(QTreeWidget* treeWidget)
	{
		QTreeWidgetItem * massaction = new QTreeWidgetItem;
		massaction->setText(0,"Mass Action");
		massaction->setToolTip(0,"evolve networks with mass action kinetics");
		treeWidget->addTopLevelItem(massaction);
		
		QCheckBox * checkBox;
		treeWidget->setItemWidget(massaction,1,checkBox = new QCheckBox);
		checkBox->setChecked(mass_action_prob > 0.0);
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
		doubleSpinBox->setValue(mass_action_prob);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setMassActionProb(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Uni-Uni");
		child->setToolTip(0,"Proportion of the reactions that will have one reactant and product");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(uni_uni);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setUniUni(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Uni-Bi");
		child->setToolTip(0,"Proportion of the reactions that will have one reactant and two products");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(uni_bi);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setUniBi(double)));		
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Bi-Uni");
		child->setToolTip(0,"Proportion of the reactions that will have two reactants and one product");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(bi_uni);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setBiUni(double)));		
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent Bi-Bi");
		child->setToolTip(0,"Proportion of the reactions that will have two reactants and products");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(bi_bi);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setBiBi(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent No reactants");
		child->setToolTip(0,"Proportion of the reactions that represent flux into the system (no reactant). This proportion is independent of the other proportions (uni-uni, bi-uni, etc.)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(no_reactant);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setNoReactant(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Percent No products");
		child->setToolTip(0,"Proportion of the reactions that represent flow out of the system (no product). This proportion is independent of the other proportions (uni-uni, bi-uni, etc.)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(no_product);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setNoProduct(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Rate constant range");
		child->setToolTip(0,"Minimum and maximum value for a reaction rate constant");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(ma_min_constant);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_min_constant(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(ma_max_constant);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_max_constant(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate rate constant");
		child->setToolTip(0,"Probability of mutating a rate constant during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(ma_mutate_constants);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_constants(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate remove reaction");
		child->setToolTip(0,"Probability of removing a reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(ma_mutate_remove_reaction);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_remove_reaction(double)));
		
		massaction->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate add reaction");
		child->setToolTip(0,"Probability of adding a new reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(ma_mutate_add_reaction);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_ma_mutate_add_reaction(double)));
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
		checkBox->setChecked(enzyme_prob > 0.0);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useEnzyme(bool)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be enzymatic networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_prob);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setEnzymeProb(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Kcat range");
		child->setToolTip(0,"Minimum and maximum for the Kcat parameter");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_min_kcat);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_kcat(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_max_kcat);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_kcat(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Keq range (log)");
		child->setToolTip(0,"Minimum and maximum for the Keq parameter in log scale (base 2)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_min_log_keq);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_log_keq(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_max_log_keq);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_log_keq(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Alpha range (log)");
		child->setToolTip(0,"Minimum and maximum for the Alpha parameter in log scale (base 2)");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_min_alpha);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_alpha(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,32.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_max_alpha);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_alpha(double)));	
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Hill range");
		child->setToolTip(0,"Minimum and maximum for the Hill coefficient");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_min_h);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_h(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_max_h);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_h(double)));		
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"S-half range");
		child->setToolTip(0,"Minimum and maximum for the substrate half saturation value");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_min_s_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_s_half(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_max_s_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_s_half(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"S-half range");
		child->setToolTip(0,"Minimum and maximum for the substrate half saturation value");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_min_p_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_min_p_half(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_max_p_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_max_p_half(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate enzyme");
		child->setToolTip(0,"Probability of swapping an enzyme during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(enzyme_mutate_enzyme);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_enzyme(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Kcat");
		child->setToolTip(0,"Probability of mutating a Kcat (catalytic constant) during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_k_cat);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_k_cat(double)));		
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Keq");
		child->setToolTip(0,"Probability of mutating the equilibrium constant during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_k_eq);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_k_eq(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Alpha");
		child->setToolTip(0,"Probability of mutating the Alpha value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_alpha);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_alpha(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Hill");
		child->setToolTip(0,"Probability of mutating the Hill coefficient during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_h);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_h(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate S-half");
		child->setToolTip(0,"Probability of mutating the substrate half-saturation point during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_s_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_s_half(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate P-half");
		child->setToolTip(0,"Probability of mutating the product half-saturation point during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_p_half);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_p_half(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate remove reaction");
		child->setToolTip(0,"Probability of removing a reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_remove);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_remove(double)));
		
		enzyme->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate add reaction");
		child->setToolTip(0,"Probability of adding a new reacton during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(enzyme_mutate_add);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_enzyme_mutate_add(double)));
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
		checkBox->setChecked(protein_net_prob > 0.0);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useProtein(bool)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be protein interaction networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(protein_net_prob);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setProteinProb(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Km range");
		child->setToolTip(0,"The minimum and maximum Michaelis-Menten constant for initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_min_ka);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_min_ka(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_max_ka);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_max_ka(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Vmax range");
		child->setToolTip(0,"The minimum and maximum Vmax constant for initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_min_vmax);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_min_vmax(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_max_vmax);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_max_vmax(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Total range");
		child->setToolTip(0,"The minimum and maximum total concentration for proteins in the initial networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_min_total);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_min_total(double)));		
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_max_total);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_max_total(double)));		
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate regulation");
		child->setToolTip(0,"Probability of rewiring the network during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_mutate_rewire);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_rewire(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate parameter");
		child->setToolTip(0,"Probability of changing a parameter during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_mutate_parameter);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_parameter(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate total");
		child->setToolTip(0,"Probability of changing a conservation law (total concentration) during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_mutate_total);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_total(double)));
		
		protein->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Add/remove proteins");
		child->setToolTip(0,"Probability of adding or removing proteins in the network during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(prot_mutate_addremove);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_prot_mutate_addremove(double)));
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
		checkBox->setChecked(grn_prob > 0.0);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(useGRN(bool)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Probability");
		child->setToolTip(0,"Proportion of the networks in the population that will be genetic networks");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);
		doubleSpinBox->setValue(grn_prob);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setGRNProb(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Resouce inflow");
		child->setToolTip(0,"Use a non-zero value here to model resource consumption. This value determines the constant inflow of resources from the environment.");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(grn_init_inflow);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_inflow(double)));		
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Resouce consumption");
		child->setToolTip(0,"Use a non-zero value here to model resource consumption. This value determines the amount of resources consumed due to protein production.");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_init_cost_per_protein);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_init_cost_per_protein(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Complex size range");
		child->setToolTip(0,"The minimum and maximum number of transcription factors that can form a complex");
		treeWidget->setItemWidget(child,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,10);
		intSpinBox->setSingleStep(1);		
		intSpinBox->setValue(grn_min_complex_size);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_grn_min_complex_size(int)));
		treeWidget->setItemWidget(child,2,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,10);
		intSpinBox->setSingleStep(1);		
		intSpinBox->setValue(grn_max_complex_size);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_grn_max_complex_size(int)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Ka range");
		child->setToolTip(0,"Minimum and maximum association constant for transcription factors");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_min_Ka);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_min_Ka(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_max_Ka);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_max_Ka(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Vmax range");
		child->setToolTip(0,"Minimum and maximum Vmax for protein production rate");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(grn_min_Vmax);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_min_Vmax(double)));
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);	
		doubleSpinBox->setValue(grn_max_Vmax);	
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_max_Vmax(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Degradation range");
		child->setToolTip(0,"Minimum and maximum degradation/dilution rate for proteins");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_min_degradation);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_min_degradation(double)));		
		treeWidget->setItemWidget(child,2,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,10000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_max_degradation);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_max_degradation(double)));		
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Ka");
		child->setToolTip(0,"Probability of changing a Ka value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_mutate_Ka);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_Ka(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Vmax");
		child->setToolTip(0,"Probability of changing a Vmax value during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_mutate_Vmax);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_Vmax(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Mutate Complex");
		child->setToolTip(0,"Probability of changing the proteins inside a complex during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_mutate_complex);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_complex(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Add genes");
		child->setToolTip(0,"Probability of adding a new gene during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_mutate_add_gene);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_add_gene(double)));
		
		grn->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Remove genes");
		child->setToolTip(0,"Probability of removing a gene during a mutation event");
		treeWidget->setItemWidget(child,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(grn_mutate_remove_gene);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_grn_mutate_remove_gene(double)));
	}
	
	
	QWidget * MainWindow::setupNetworkOptions()
	{
		QTreeWidget * treeWidget = new QTreeWidget();
		treeWidget->setColumnCount(3);
		treeWidget->setColumnWidth(0,200);
		treeWidget->setHeaderLabels( QStringList() << "network/property" << "value/min" << "max");
		
		setupMassActionNetwork(treeWidget);
		setupEnzymeNetwork(treeWidget);
		setupProteinNetwork(treeWidget);
		setupGeneticNetwork(treeWidget);
		
		
		QCheckBox * checkBox;
		QDoubleSpinBox * doubleSpinBox;
		QSpinBox * intSpinBox;
		QTreeWidgetItem * option;
		
		option = new QTreeWidgetItem;
		option->setText(0,"Num. Molecular Species");
		option->setToolTip(0,"Minimum and maximum number of molecular species is the evolved networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);	
		intSpinBox->setValue(species_min);	
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_min_species(int)));
		treeWidget->setItemWidget(option,2,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);	
		intSpinBox->setValue(species_max);	
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_max_species(int)));
		
		option = new QTreeWidgetItem;
		option->setText(0,"Num. Reactions");
		option->setToolTip(0,"Minimum and maximum number of reactions is the evolved networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);		
		intSpinBox->setValue(reactions_min);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_min_reactions(int)));
		treeWidget->setItemWidget(option,2,intSpinBox = new QSpinBox);
		intSpinBox->setRange(0,100000);
		intSpinBox->setSingleStep(1);		
		intSpinBox->setValue(reactions_max);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(set_max_reactions(int)));
		
		option = new QTreeWidgetItem;
		option->setText(0,"Crossover Rate");
		option->setToolTip(0,"Probability of a crossover event when generating the next generation of networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(crossover_rate);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_crossover_rate(double)));
		
		option = new QTreeWidgetItem;
		option->setText(0,"Avg. initial values");
		option->setToolTip(0,"Average concentration values of molecules in the initial networks");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,100000.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(init_iv);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_init_iv(double)));
		
		option = new QTreeWidgetItem;
		option->setText(0,"Mutate initial values");
		option->setToolTip(0,"Probability of mutating the concentration values during a mutation event");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,doubleSpinBox = new QDoubleSpinBox);
		doubleSpinBox->setRange(0.0,1.0);
		doubleSpinBox->setDecimals(3);
		doubleSpinBox->setSingleStep(0.01);		
		doubleSpinBox->setValue(mutate_iv);		
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(set_mutate_iv(double)));
		
		option = new QTreeWidgetItem;
		option->setText(0,"Lineage tracking");
		option->setToolTip(0,"Track the original parent(s) for each network, and record this table in the log file");
		treeWidget->addTopLevelItem(option);
		treeWidget->setItemWidget(option,1,checkBox = new QCheckBox);
		checkBox->setChecked(lineageTracking);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(set_lineageTracking(bool)));
		
		
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
		treeWidget->setItemWidget(code,1,lineEdit = new QLineEdit(codeFile));
		connect(lineEdit,SIGNAL(textEdited (const QString &)),this,SLOT(setCodeFile(const QString&)));
		
		QTreeWidgetItem * log = new QTreeWidgetItem;
		treeWidget->addTopLevelItem(log);
		treeWidget->setItemWidget(log,1,lineEdit = new QLineEdit(logFile));
		connect(lineEdit,SIGNAL(textEdited (const QString &)),this,SLOT(setLogFile(const QString&)));
		
		QTreeWidgetItem * runs = new QTreeWidgetItem;
		runs->setText(0,"Runs");
		runs->setToolTip(0,"The number of times to repeat the evolution experiment. Each result will be different, unless the same seed is used");
		treeWidget->addTopLevelItem(runs);
		treeWidget->setItemWidget(runs,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		intSpinBox->setValue(this->runs);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setRuns(int)));
		
		QTreeWidgetItem * popSz = new QTreeWidgetItem;
		popSz->setText(0,"Population Size");
		popSz->setToolTip(0,"The final population size for each run");
		treeWidget->addTopLevelItem(popSz);
		treeWidget->setItemWidget(popSz,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		intSpinBox->setValue(this->popSz);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setPopSz(int)));
		
		QTreeWidgetItem * initPopSz = new QTreeWidgetItem;
		initPopSz->setText(0,"Initial Population Size");
		initPopSz->setToolTip(0,"The initial parent population size for each run");
		treeWidget->addTopLevelItem(initPopSz);
		treeWidget->setItemWidget(initPopSz,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		intSpinBox->setValue(this->initPopSz);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setInitPopSz(int)));
		
		QTreeWidgetItem * generation = new QTreeWidgetItem;
		generation->setText(0,"Generations");
		generation->setToolTip(0,"The number of iterations to run the genetic algorithm during each evolution experiment");
		treeWidget->addTopLevelItem(generation);
		treeWidget->setItemWidget(generation,1,intSpinBox = new QSpinBox);
		intSpinBox->setRange(1,1000000);
		intSpinBox->setSingleStep(1);
		connect(intSpinBox,SIGNAL(valueChanged(int)),this,SLOT(setGenerations(int)));
		intSpinBox->setValue(generations);
		
		QTreeWidgetItem * stopcrit = new QTreeWidgetItem;
		stopcrit->setText(0,"Stop criterion");
		stopcrit->setToolTip(0,"Stop the evolution when best individual reaches this fitness level. Use 0 to disable.");
		treeWidget->addTopLevelItem(stopcrit);
		treeWidget->setItemWidget(stopcrit,1,doubleSpinBox = new QDoubleSpinBox);
		//doubleSpinBox->setRange(1,1000000);
		doubleSpinBox->setSingleStep(0.01);
		doubleSpinBox->setDecimals(5);
		connect(doubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(setMaxFitness(double)));
		doubleSpinBox->setValue(max_fitness);
		
		log->setText(0,"Log file");
		log->setToolTip(0,"Keep a log file with information about each generation and the final results");
		
		QTreeWidgetItem * log1, * log2, *child;
		log->addChild(log1 = new QTreeWidgetItem);
		log1->setText(0,"Track information during each generation");
		log1->setToolTip(0,"Check the information you would like to record in the log file during each generation of the genetic algorithm.");
		
		log->addChild(log2 = new QTreeWidgetItem);
		log2->setText(0,"Report final results");
		log2->setToolTip(0,"Check the information you would like to record at the end of each run, where one run is one evolution experiment");
		
		QCheckBox * checkBox;
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best fitness score");
		child->setToolTip(0,"Report the fitness score of the best network in this generation.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkFitness1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkFitness1(bool)));
		
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network script");
		child->setToolTip(0,"Print the script for the best network during each generation. Warning: expect a bloated log file");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkScript1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkScript1(bool)));
		
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network's size");
		child->setToolTip(0,"Record the best network's size during each generation.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkSize1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkSize1(bool)));
		
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network's parents");
		child->setToolTip(0,"Record the best network's lineage during each generation.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkLineage1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkLineage1(bool)));
		
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"All fitness scores");
		child->setToolTip(0,"Record all the networks' fitness scores during each generation.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(allFitness1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setAllFitness1(bool)));
		
		log1->addChild(child = new QTreeWidgetItem);
		child->setText(0,"All parents");
		child->setToolTip(0,"Record all parents. Enable this option if you are interested in the evolutionary dynamics of the population.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(allNetworkLineage1);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setAllNetworkLineage1(bool)));
		
		//---
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best fitness score");
		child->setToolTip(0,"Report the fitness score of the best network at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkFitness2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkFitness2(bool)));
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network script");
		child->setToolTip(0,"Print the script for the best network at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkScript2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkScript2(bool)));
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network's size");
		child->setToolTip(0,"Record the best network's size  at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkSize2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkSize2(bool)));
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Best network's parents");
		child->setToolTip(0,"Record the best network's lineage at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(bestNetworkLineage2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setBestNetworkLineage2(bool)));
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"All fitness scores");
		child->setToolTip(0,"Record all the networks' fitness scores at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(allFitness2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setAllFitness2(bool)));
		
		log2->addChild(child = new QTreeWidgetItem);
		child->setText(0,"All parents");
		child->setToolTip(0,"Record all parents at the end of each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(allNetworkLineage2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(setAllNetworkLineage2(bool)));
		
		log->addChild(child = new QTreeWidgetItem);
		child->setText(0,"Save seeds");
		child->setToolTip(0,"Record the random number generator's seeds that were used for each run.");
		treeWidget->setItemWidget(child,1,checkBox = new QCheckBox);
		checkBox->setChecked(seeds2);
		connect(checkBox,SIGNAL(toggled(bool)),this,SLOT(showSeed(bool)));
		
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
		fitnessComboBox = comboBox();
		connect(fitnessComboBox,SIGNAL(activated(const QString&)),this,SLOT(fitnessSelected(const QString&)));
		
		codeEditor = new Tinkercell::CodeEditor;
		codeEditor->setStyleSheet("background-color: #F9F7E4");
		codeEditor->setTabStopWidth(20);
		
		CSyntaxHighlighter * syntaxHighlighter = new CSyntaxHighlighter(codeEditor->document());
		
		QVBoxLayout * layout = new QVBoxLayout;
		layout->addWidget(fitnessComboBox);
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
				tr("#include \"reactionNetwork.h\"\n\n")
				+ tr("void init()\n")
				+ tr("{\n")
				+ tr("    setNetworkTypeProbability(0,") + QString::number(mass_action_prob) + tr(");\n")
				+ tr("    setNetworkTypeProbability(1,") + QString::number(enzyme_prob) + tr(");\n")
				+ tr("    setNetworkTypeProbability(2,") + QString::number(protein_net_prob) + tr(");\n")
				+ tr("    setNetworkTypeProbability(3,") + QString::number(grn_prob) + tr(");\n\n")
				+ tr("    setDistributionOfMassActionNetwork(")
				+ QString::number(uni_uni) + tr(",")
				+ QString::number(uni_bi) + tr(",")
				+ QString::number(bi_uni) + tr(",")
				+ QString::number(bi_bi) + tr(",")
				+ QString::number(no_reactant) + tr(",")
				+ QString::number(no_product)
				+ tr(");\n")
				+ tr("    setRateConstantForMassActionNetwork(")
				+ QString::number(ma_min_constant) + tr(",")
				+ QString::number(ma_max_constant)
				+ tr(");\n")
				+ tr("    setMutationRatesForMassActionNetwork(")
				+ QString::number(ma_mutate_constants) + tr(",")
				+ QString::number(ma_mutate_remove_reaction) + tr(",")
				+ QString::number(ma_mutate_add_reaction)
				+ tr(");\n")
				+ tr("    setRateConstantsForProteinInteractionNetwork(")
				+ QString::number(prot_min_ka) + tr(",")
				+ QString::number(prot_max_ka) + tr(",")
				+ QString::number(prot_min_vmax) + tr(",")
				+ QString::number(prot_max_vmax) + tr(",")
				+ QString::number(prot_min_total) + tr(",")
				+ QString::number(prot_max_total)
				+ tr(");\n")
				+ tr("    setMutationRatesForProteinInteractionNetwork(")
				+ QString::number(prot_mutate_rewire) + tr(",")
				+ QString::number(prot_mutate_parameter) + tr(",")
				+ QString::number(prot_mutate_total) + tr(",")
				+ QString::number(prot_mutate_addremove)
				+ tr(");\n")
				+ tr("    setRateConstantsForEnzymeNetwork(")
				+ QString::number(enzyme_min_kcat) + tr(",")
				+ QString::number(enzyme_max_kcat) + tr(",")
				+ QString::number(enzyme_min_log_keq) + tr(",")
				+ QString::number(enzyme_max_log_keq) + tr(",")
				+ QString::number(enzyme_min_alpha) + tr(",")
				+ QString::number(enzyme_max_alpha) + tr(",")
				+ QString::number(enzyme_min_h) + tr(",")
				+ QString::number(enzyme_max_h) + tr(",")
				+ QString::number(enzyme_min_s_half) + tr(",")
				+ QString::number(enzyme_max_s_half) + tr(",")
				+ QString::number(enzyme_min_p_half) + tr(",")
				+ QString::number(enzyme_max_p_half)
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
				+ QString::number(grn_min_complex_size) + tr(",")
				+ QString::number(grn_max_complex_size) + tr(",")
				+ QString::number(grn_min_Ka) + tr(",")
				+ QString::number(grn_max_Ka) + tr(",")
				+ QString::number(grn_min_Vmax) + tr(",")
				+ QString::number(grn_max_Vmax) + tr(",")
				+ QString::number(grn_min_degradation) + tr(",")
				+ QString::number(grn_max_degradation)
				+ tr(");\n")
				+ tr("    setMutationRatesForGeneRegulationNetwork(")
				+ QString::number(grn_mutate_Ka) + tr(",")
				+ QString::number(grn_mutate_Vmax) + tr(",")
				+ QString::number(grn_mutate_complex) + tr(",")
				+ QString::number(grn_mutate_add_gene) + tr(",")
				+ QString::number(grn_mutate_remove_gene)
				+ tr(");\n")
				+ tr("    setNetworkSize(")
				+ QString::number(species_min) + tr(",")
				+ QString::number(species_max) + tr(",")
				+ QString::number(reactions_min) + tr(",")
				+ QString::number(reactions_max) 
				+ tr(");\n")
				+ tr("    setCrossoverRate(")
				+ QString::number(crossover_rate) 
				+ tr(");\n")
				+ tr("    setAverageInitialValue(")
				+ QString::number(init_iv) 
				+ tr(");\n")
				+ tr("    setMutationRateOfInitialValues(")
				+ QString::number(mutate_iv) 
				+ tr(");\n")
				+ tr("    configureContinuousLog(") 
				+ QString::number(bestNetworkFitness1) + tr(",")
				+ QString::number(bestNetworkScript1) + tr(",")
				+ QString::number(bestNetworkSize1) + tr(",")
				+ QString::number(bestNetworkLineage1) + tr(",")
				+ QString::number(allFitness1) + tr(",")
				+ QString::number(allNetworkLineage1) + tr(");\n")
				+ tr("    configureFinalLog(")
				+ QString::number(bestNetworkFitness2) + tr(",")
				+ QString::number(bestNetworkScript2) + tr(",")
				+ QString::number(bestNetworkSize2) + tr(",")
				+ QString::number(bestNetworkLineage2) + tr(",")
				+ QString::number(allFitness2) + tr(",")
				+ QString::number(allNetworkLineage2) + tr(",")
				+ QString::number(seeds2) + tr(");\n");
		
		if (logFile.isEmpty())
			s += tr("    disableLogFile();\n");
		else
			s += tr("    enableLogFile(\"") + logFile + tr("\");\n");
		
		if (lineageTracking)
			s += tr("    lineageTrackingON();\n");
		else
			s += tr("    lineageTrackingOFF();\n");
			
		if (!seeds.isEmpty() && seeds.split(tr(",")).size() == 4)
		{
			s +=  tr("    unsigned long long a[] = {") + seeds + tr("};\n")
				+ tr("    void setMTseeds(") + s + tr(");\n");
		}
		
		s += tr("}\n");
		
		return s;
	}
	
	QString MainWindow::callbackFunction()
	{
		if (max_fitness > 0)
		{
			return 
			tr("int callback(int iter, GApopulation P, int popSz)\n")
			+ tr("{\n")
			+ tr("    return (fitness(P[0]) > ")
			+ QString::number(max_fitness)
			+ tr(");\n}\n");
		}
		
		return 
			tr("int callback(int iter, GApopulation P, int popSz)\n")
			+ tr("{\n")
			+ tr("    return 0;\n}\n");
	}
	
	QString MainWindow::mainFunction()
	{
		QString s = tr("int main()\n")
			+ tr("{\n")
			+ tr("    int runs = ") + QString::number(runs) + tr(";\n")
			+ tr("    int generations = ") + QString::number(generations) + tr(";\n")
			+ tr("    int popSz = ") + QString::number(popSz) + tr(";\n")
			+ tr("    int initPopSz = ") + QString::number(initPopSz) + tr(";\n")
			+ tr("    GApopulation P;\n")
			+ tr("    int i,j;\n\n")
			+ tr("    setFitnessFunction(&fitness);\n")
			+ tr("    init();\n")
			+ tr("    for (i=0; i < runs; ++i)\n")
			+ tr("    {\n")
			+ tr("       printf(\"run #%i\",i+1);\n")
			+ tr("       P = evolveNetworks(initPopSz,popSz,generations,&callback);\n")
			+ tr("       GAfree(P);\n")
			+ tr("    }\n")
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
			callbackFunction() + tr("\n");

		qfile.close();
		
		proc.setWorkingDirectory(appDir);
		
#ifdef Q_WS_WIN
		
		proc.start(compileCommand);
		proc.waitForFinished();
		proc.start(
					tr("RunEvolution.exe temp.dll ")
					+ tr(" ") + QString::number(generations) 
					+ tr(" ") + QString::number(initPopSz) 
					+ tr(" ") + QString::number(popSz) );
		//proc.waitForFinished();
		
#else
#ifdef Q_WS_MAC

		proc.start(compileCommand);
		proc.waitForFinished();
		proc.start(
					tr("./RunEvolution temp.dylib")
					+ tr(" ") + QString::number(generations) 
					+ tr(" ") + QString::number(initPopSz) 
					+ tr(" ") + QString::number(popSz) );
		//proc.waitForFinished();
		
#else
		proc.start(compileCommand);
		proc.waitForFinished();
		proc.start(
					tr("./RunEvolution temp.so ")
					+ tr(" ") + QString::number(generations) 
					+ tr(" ") + QString::number(initPopSz) 
					+ tr(" ") + QString::number(popSz) );
		//proc.waitForFinished();
		
#endif
#endif

	}
	
	void MainWindow::reset()
	{
		/***************************/
		/*****initial values*******/
		/**************************/
		
		mass_action_prob = 0.0;
		enzyme_prob = 0.0;
		protein_net_prob = 0.0;
		grn_prob = 1.0;
		
		uni_uni = uni_bi = bi_uni = bi_bi = no_reactant = no_product = 0.2;
		ma_min_constant = 0.01;
		ma_max_constant = 100.0;
		ma_mutate_constants = 0.5;
		ma_mutate_remove_reaction = ma_mutate_add_reaction = 0.25;
		
		enzyme_max_kcat = 100.0;
		enzyme_max_log_keq = enzyme_max_h = 4.0;
		enzyme_max_alpha = 4.0;
		enzyme_max_s_half = enzyme_max_p_half = 20.0;
		
		enzyme_min_kcat = 0.001;
		enzyme_min_log_keq = -4.0;
		enzyme_min_h = 1.0;
		enzyme_min_alpha = -4.0;
		enzyme_min_s_half = enzyme_min_p_half = 0.001;
		
		enzyme_mutate_enzyme = enzyme_mutate_k_eq = enzyme_mutate_k_cat = enzyme_mutate_alpha = 0.1;
		enzyme_mutate_h = enzyme_mutate_s_half = enzyme_mutate_p_half = 0.1;
		enzyme_mutate_remove = enzyme_mutate_add = 0.1;		
		
		prot_min_ka = 0.001;
		prot_max_ka = 100.0;
		prot_min_vmax = 0.1;
		prot_max_vmax = 100.0;
		prot_min_total = 0.1;
		prot_max_total = 100.0;
		prot_mutate_rewire = prot_mutate_total = prot_mutate_addremove = 0.2;
		prot_mutate_parameter = 0.4;
		
		grn_init_inflow = grn_init_cost_per_protein = 0.0;
		grn_min_complex_size = 1;
		grn_max_complex_size = 4;
		grn_min_Ka = 0.001;
		grn_max_Ka = 100.0;
		grn_min_Vmax = 0.1;
		grn_max_Vmax = 20.0;
		grn_min_degradation = 0.1;
		grn_max_degradation = 10.0;
		grn_mutate_Ka = grn_mutate_Vmax = grn_mutate_complex = grn_mutate_add_gene = grn_mutate_remove_gene = 0.2;
		
		codeFile = tr("run_net_ga.c");
		logFile = tr("evolution.log");
		runs = 10;
		popSz = 200;
		initPopSz = 200;
		generations = 50;
		
		species_min = 4;
		species_max = 12;
		reactions_min = 5;
		reactions_max = 24;
		crossover_rate = init_iv = 1.0;
		mutate_iv = 0.05;
		lineageTracking = true;
		
		bestNetworkFitness1 = 1;
		bestNetworkScript1 = 0;
		bestNetworkSize1 = 1;
		bestNetworkLineage1 = 0;
		allFitness1 = 0;
		allNetworkLineage1 = 1;
			
		bestNetworkFitness2 = 1;
		bestNetworkScript2 = 1;
		bestNetworkSize2 = 1;
		bestNetworkLineage2 = 0;
		allFitness2 = 1;
		allNetworkLineage2 = 1;
		seeds2 = 1;
		
		max_fitness = 0.0;
		
		/***********************/
		
		QSettings settings("UWashington","NetworkEvolutionLib");
		settings.beginGroup("parameters");
		
		uni_uni = settings.value("uni_uni",uni_uni).toDouble();
		uni_bi = settings.value("uni_bi",uni_bi).toDouble();
		bi_uni = settings.value("bi_uni",bi_uni).toDouble();
		bi_uni = settings.value("bi_bi",bi_uni).toDouble();
		no_reactant = settings.value("no_reactant",no_reactant).toDouble();
		no_product = settings.value("no_product",no_product).toDouble();
		ma_min_constant = settings.value("ma_min_constant",ma_min_constant).toDouble();
		ma_max_constant = settings.value("ma_max_constant",ma_max_constant).toDouble();
		ma_mutate_constants = settings.value("ma_mutate_constants",ma_mutate_constants).toDouble();
		ma_mutate_remove_reaction = settings.value("ma_mutate_remove_reaction",ma_mutate_remove_reaction).toDouble();
		ma_mutate_add_reaction = settings.value("ma_mutate_add_reaction",ma_mutate_add_reaction).toDouble();
		
		prot_min_ka = settings.value("prot_min_ka",prot_min_ka).toDouble();
		prot_min_vmax = settings.value("prot_min_vmax",prot_min_vmax).toDouble();
		prot_min_total = settings.value("prot_min_total",prot_min_total).toDouble();
		prot_max_ka = settings.value("prot_max_ka",prot_max_ka).toDouble();
		prot_max_vmax = settings.value("prot_max_vmax",prot_max_vmax).toDouble();
		prot_max_total = settings.value("prot_max_total",prot_max_total).toDouble();
		
		prot_mutate_rewire = settings.value("prot_mutate_rewire",prot_mutate_rewire).toDouble();
		prot_mutate_parameter = settings.value("prot_mutate_parameter",prot_mutate_parameter).toDouble();
		prot_mutate_total = settings.value("prot_mutate_total",prot_mutate_total).toDouble();
		prot_mutate_addremove = settings.value("prot_mutate_addremove",prot_mutate_addremove).toDouble();
		
		enzyme_min_kcat = settings.value("enzyme_min_kcat",enzyme_min_kcat).toDouble();
		enzyme_min_log_keq = settings.value("enzyme_min_log_keq",enzyme_min_log_keq).toDouble();
		enzyme_min_alpha = settings.value("enzyme_min_alpha",enzyme_min_alpha).toDouble();
		enzyme_min_h = settings.value("enzyme_min_h",enzyme_min_h).toDouble();
		enzyme_min_s_half = settings.value("enzyme_min_s_half",enzyme_min_s_half).toDouble();
		enzyme_min_p_half = settings.value("enzyme_min_p_half",enzyme_min_p_half).toDouble();
		enzyme_max_kcat = settings.value("enzyme_max_kcat",enzyme_max_kcat).toDouble();
		enzyme_max_log_keq = settings.value("enzyme_max_log_keq",enzyme_max_log_keq).toDouble();
		enzyme_max_alpha = settings.value("enzyme_max_alpha",enzyme_max_alpha).toDouble();
		enzyme_max_h = settings.value("enzyme_max_h",enzyme_max_h).toDouble();
		enzyme_max_s_half = settings.value("enzyme_max_s_half",enzyme_max_s_half).toDouble();
		enzyme_max_p_half = settings.value("enzyme_max_p_half",enzyme_max_p_half).toDouble();
		
		enzyme_mutate_enzyme = settings.value("enzyme_mutate_enzyme",enzyme_mutate_enzyme).toDouble();
		enzyme_mutate_k_cat = settings.value("enzyme_mutate_k_cat",enzyme_mutate_k_cat).toDouble();
		enzyme_mutate_k_cat = settings.value("enzyme_mutate_k_eq",enzyme_mutate_k_eq).toDouble();
		enzyme_mutate_h = settings.value("enzyme_mutate_h",enzyme_mutate_h).toDouble();
		enzyme_mutate_s_half = settings.value("enzyme_mutate_s_half",enzyme_mutate_s_half).toDouble();
		enzyme_mutate_p_half = settings.value("enzyme_mutate_p_half",enzyme_mutate_p_half).toDouble();
		enzyme_mutate_remove = settings.value("enzyme_mutate_remove",enzyme_mutate_remove).toDouble();
		enzyme_mutate_add = settings.value("enzyme_mutate_add",enzyme_mutate_add).toDouble();
		
		grn_init_inflow = settings.value("grn_init_inflow",grn_init_inflow).toDouble();
		grn_init_cost_per_protein = settings.value("grn_init_cost_per_protein",grn_init_cost_per_protein).toDouble();
		grn_min_complex_size = settings.value("grn_min_complex_size",grn_min_complex_size).toInt();
		grn_min_Ka = settings.value("grn_min_Ka",grn_min_Ka).toDouble();
		grn_min_Vmax = settings.value("grn_min_Vmax",grn_min_Vmax).toDouble();
		grn_min_degradation = settings.value("grn_min_degradation",grn_min_degradation).toDouble();
		grn_max_complex_size = settings.value("grn_max_complex_size",grn_max_complex_size).toInt();
		grn_max_Ka = settings.value("grn_max_Ka",grn_max_Ka).toDouble();
		grn_max_Vmax = settings.value("grn_max_Vmax",grn_max_Vmax).toDouble();
		grn_max_degradation = settings.value("grn_max_degradation",grn_max_degradation).toDouble();
		
		grn_mutate_Ka = settings.value("grn_mutate_Ka",grn_mutate_Ka).toDouble();
		grn_mutate_Vmax = settings.value("grn_mutate_Vmax",grn_mutate_Vmax).toDouble();
		grn_mutate_complex = settings.value("grn_mutate_complex",grn_mutate_complex).toDouble();
		grn_mutate_add_gene = settings.value("grn_mutate_add_gene",grn_mutate_add_gene).toDouble();
		grn_mutate_remove_gene = settings.value("grn_mutate_remove_gene",grn_mutate_remove_gene).toDouble();
		
		species_min = settings.value("species_min",species_min).toInt();
		species_max = settings.value("species_max",species_max).toInt();
		reactions_min = settings.value("reactions_min",reactions_min).toInt();
		reactions_max = settings.value("reactions_max",reactions_max).toInt();
		crossover_rate = settings.value("crossover_rate",crossover_rate).toDouble();
		init_iv = settings.value("init_iv",init_iv).toDouble();
		mutate_iv = settings.value("mutate_iv",mutate_iv).toDouble();
		lineageTracking = settings.value("lineageTracking",lineageTracking).toBool();
		
		mass_action_prob = settings.value("mass_action_prob",mass_action_prob).toDouble();
		enzyme_prob = settings.value("enzyme_prob",enzyme_prob).toDouble();
		protein_net_prob = settings.value("protein_net_prob",protein_net_prob).toDouble();
		grn_prob = settings.value("grn_prob",grn_prob).toDouble();
		
		codeFile = settings.value("codeFile",codeFile).toString();
		logFile = settings.value("logFile",logFile).toString();
		runs = settings.value("runs",runs).toInt();
		generations = settings.value("generations",generations).toInt();
		popSz = settings.value("popSz",popSz).toInt();
		initPopSz = settings.value("initPopSz",initPopSz).toInt();
		
		bestNetworkFitness1 = settings.value("bestNetworkFitness1",bestNetworkFitness1).toInt();
		bestNetworkScript1 = settings.value("bestNetworkScript1",bestNetworkScript1).toInt();
		bestNetworkSize1 = settings.value("bestNetworkSize1",bestNetworkSize1).toInt();
		bestNetworkLineage1 = settings.value("bestNetworkLineage1",bestNetworkLineage1).toInt();
		allFitness1 = settings.value("allFitness1",allFitness1).toInt();
		allNetworkLineage1 = settings.value("allNetworkLineage1",allNetworkLineage1).toInt();
		
		bestNetworkFitness2 = settings.value("bestNetworkFitness2",bestNetworkFitness2).toInt();
		bestNetworkScript2 = settings.value("bestNetworkScript2",bestNetworkScript2).toInt();
		bestNetworkSize2 = settings.value("bestNetworkSize2",bestNetworkSize2).toInt();
		bestNetworkLineage2 = settings.value("bestNetworkLineage2",bestNetworkLineage2).toInt();
		allFitness2 = settings.value("allFitness2",allFitness2).toInt();
		allNetworkLineage2 = settings.value("allNetworkLineage2",allNetworkLineage2).toInt();
		seeds2 = settings.value("seeds2",seeds2).toInt();
		
		max_fitness = settings.value("max_fitness",max_fitness).toDouble();
		
#ifdef Q_WS_WIN
		compileCommand = tr("gcc --shared -o temp.dll -I..\\..\\GA -I..\\..\\simulation temp.c -L./ -lnetga");
#else
#ifdef Q_WS_MAC
		compileCommand = tr("gcc --shared -o temp.dylib -I../../GA -I../../simulation temp.c -L./ -lnetga");
#else
		compileCommand = tr("gcc --shared -o temp.so -I../../GA -I../../simulation temp.c -L./ -lnetga");
#endif
#endif
		
		settings.endGroup();
	}
	
	
	QComboBox* MainWindow::comboBox()
	{
		QComboBox * comboBox = new QComboBox;
		comboBox->addItem("NEW");
		
		QString appDir = QCoreApplication::applicationDirPath();
		
		QDir dir(appDir + tr("/FitnessFunctions/"));
		dir.setNameFilters (QStringList() << "*.c");
		dir.setSorting(QDir::Name);

		QFileInfoList list = dir.entryInfoList();
		 
		for (int i = 0; i < list.size(); ++i) 
		{
			QFileInfo fileInfo = list.at(i);
			comboBox->addItem(fileInfo.baseName());
		}
		
		return comboBox;
	}
	
	void MainWindow::fitnessSelected(const QString& s)
	{
		if (fitnessComboBox->currentIndex() == 0)
		{
			clear();
			return;
		}
			
		QString appDir = QCoreApplication::applicationDirPath();
		QFile file(appDir + tr("/FitnessFunctions/") + s + tr(".c"));
	
		if (file.open(QFile::ReadOnly | QFile::Text))
		{
			codeEditor->setPlainText(QString(file.readAll()));
			file.close();
		}
	}
	
	void MainWindow::clear()
	{
		QString s = tr("double fitness(GAindividual x)\n{\n    return mtrand();\n }\n");
		codeEditor->setPlainText(s);
	}
	
	void MainWindow::save()
	{
		QString appDir = QCoreApplication::applicationDirPath();
		
		QString fileName =
					QFileDialog::getSaveFileName(this, tr("Save Current Model"),
						appDir + tr("/FitnessFunctions/"),
						tr("C Files (*.c)"));
		if (fileName.isEmpty())
			return;
		
		QFile file(fileName);
		if (file.open(QFile::WriteOnly | QFile::Text))
		{
			fitnessComboBox->addItem(QFileInfo(file).baseName());
			file.write( codeEditor->toPlainText().toAscii() );
			file.close();
		}
	}
	
}
