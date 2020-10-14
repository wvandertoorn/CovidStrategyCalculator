#include "mainwindow.h"

#include <QGridLayout>
#include <QGroupBox>
#include <QPushButton>
#include <QDesktopWidget>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    QSize size = QDesktopWidget().availableGeometry(this).size();
    resize(size);

    tab = new QTabWidget;
    tab->addTab(initialize_tab_input(), tr("Input"));
    tab->addTab(initialize_tab_parameters(), tr("Parameters"));
    tab->setCurrentIndex(0);

    chart = new QWidget;
    log = new QWidget;

    QScrollArea* scroll_tab = new QScrollArea(this);
    scroll_tab->setWidget(tab);
    scroll_tab->setWidgetResizable(true);
    scroll_tab->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Maximum);

    QScrollArea* scroll_log = new QScrollArea(this);
    scroll_log->setWidget(log);
    scroll_log->setWidgetResizable(true);
    scroll_log->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

    QScrollArea* scroll_chart = new QScrollArea(this);
    scroll_chart->setWidget(chart);
    scroll_chart->setWidgetResizable(true);
    scroll_chart->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->addWidget(scroll_tab, 0,0,1,1);
    mainLayout->addWidget(scroll_log, 0,1,1,1);
    mainLayout->addWidget(scroll_chart, 1,0,2,2);

    QWidget *centralWidget = new QWidget;
    centralWidget->setLayout(mainLayout);

    this->setCentralWidget(centralWidget);
}

QDoubleSpinBox* MainWindow::create_parameter_DoubleSpinBox(QWidget *parent, double min, double max, int dec, double val)
{
    QDoubleSpinBox *box = new QDoubleSpinBox(parent);
    box->setMinimum(min);
    box->setMaximum(max);
    box->setDecimals(dec);
    box->setValue(val);
    box->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    return box;
}

QWidget *MainWindow::initialize_tab_input()
{
    QLabel *label_mode = new QLabel(tr("Simulation start:"));
    label_mode->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    this->mode_ComboBox = new QComboBox(this);
    this->mode_ComboBox->addItems(QStringList{"exposure", "symptom onset"});
    this->mode_ComboBox->setCurrentIndex(0);

    QLabel *label_time_passed = new QLabel( ( std::string{"Time passed since "} +
                                              this->mode_ComboBox->currentText().toStdString() +
                                              std::string{" [days]:"}).c_str() );
    label_time_passed->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    this->time_passed = create_parameter_DoubleSpinBox(this, 0, 21, 0, 3);

    QLabel *label_quarantine = new QLabel(tr("Duration of quarantine [days]:"));
    label_quarantine->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    this->quarantine = create_parameter_DoubleSpinBox(this, 0, 35, 0, 10);

    this->test_days_box = new QGroupBox(this);
    this->test_days_box->setTitle(tr("Days to test on:"));
    QScrollArea* scrollArea = new QScrollArea(this);
    scrollArea->setWidget(test_days_box);
    scrollArea->setWidgetResizable(true);
    scrollArea->setHorizontalScrollBarPolicy( Qt::ScrollBarAsNeeded);
    scrollArea->setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff);
    scrollArea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);

    this->run_PushButton = new QPushButton(this);
    run_PushButton->setText(tr("Run"));
    run_PushButton->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    label_result_perc = new QLabel(tr(""));
    result_perc = new QLabel(tr(""));
    label_result_factor = new QLabel(tr(""));
    result_factor = new QLabel(tr(""));

    QGridLayout *gridLayout = new QGridLayout;
    gridLayout->addWidget(label_mode, 0, 0);
    gridLayout->addWidget(mode_ComboBox, 0, 1, Qt::AlignLeft);
    gridLayout->addWidget(label_time_passed, 1, 0);
    gridLayout->addWidget(time_passed, 1, 1, Qt::AlignLeft);
    gridLayout->addWidget(label_quarantine, 2, 0);
    gridLayout->addWidget(quarantine, 2, 1, Qt::AlignLeft);
    gridLayout->setHorizontalSpacing(10);
    gridLayout->setSizeConstraint(QLayout::SetFixedSize);

    QVBoxLayout *input_tab_layout = new QVBoxLayout;
    input_tab_layout->addItem(gridLayout);
    input_tab_layout->addWidget(scrollArea);
    input_tab_layout->addWidget(run_PushButton);
    input_tab_layout->setAlignment(Qt::AlignTop);

    connect(mode_ComboBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(mode_ComboBox_currentIndexChanged(int)) );
    connect(time_passed, SIGNAL(valueChanged(double)), this, SLOT(time_passed_valueChanged()));
    connect(quarantine, SIGNAL(valueChanged(double)), this, SLOT(quarantine_valueChanged()));
    connect(run_PushButton, SIGNAL(clicked()), this, SLOT(run_PushButton_clicked()) );

    initialize_test_date_checkboxes();

    QWidget *widget = new QWidget;
    widget->setLayout(input_tab_layout);
    return widget;
}

QWidget* MainWindow::initialize_tab_parameters()
{
    // variables
    QLabel *label_tau_inc = new QLabel(tr("Duration of incubation period [day]:"));
    QLabel *label_percentage_predetect = new QLabel(tr("Percentage thereof pre-detectable [%]:"));
    QLabel *label_tau_symp = new QLabel(tr("Duration of symptomatic period [day]:"));
    QLabel *label_asymp = new QLabel(tr("Percentage of asymptomatic infections [%]:"));
    QLabel *label_tau_post = new QLabel(tr("Duration of post-symptomatic, detectable period [day]:"));

    QLabel *label_pcr_sens = new QLabel(tr("PCR-test sensitivity [%]"));
    QLabel *label_pcr_spec = new QLabel(tr("PCR-test specificity [%]"));

    QLabel *label_mean = new QLabel(tr("Mean"));
    QLabel *label_lev = new QLabel(tr("Lower extreme value"));
    QLabel *label_uev = new QLabel(tr("Upper extreme value"));

    this->inc_lev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["inc_lev"]);
    this->inc_mean = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["inc_mean"]);
    this->inc_uev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["inc_uev"]);

    this->percentage_predetection = create_parameter_DoubleSpinBox(this, 0.01, 100, 2,
                                                                   this->default_values["percentage_predetection"]);

    this->symp_lev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["symp_lev"]);
    this->symp_mean = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["symp_mean"]);
    this->symp_uev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["symp_uev"]);

    this->percentage_asymptomatic = create_parameter_DoubleSpinBox(this, 0.01, 100, 2,
                                                                   this->default_values["percentage_asymptomatic"]);

    this->post_lev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["post_lev"]);
    this->post_mean = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["post_mean"]);
    this->post_uev = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["post_uev"]);

    this->pcr_sens = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["pcr_sens"]);
    this->pcr_spec = create_parameter_DoubleSpinBox(this, 0.01, 100, 2, this->default_values["pcr_spec"]);

    reset_PushButton = new QPushButton(this);
    reset_PushButton->setText(tr("Reset defaults"));
    reset_PushButton->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);
    connect(reset_PushButton, SIGNAL(clicked()), this, SLOT(reset_PushButton_clicked()) );

    // layout
    QGridLayout *param_tab_layout = new QGridLayout;
    param_tab_layout->setAlignment(Qt::AlignCenter);

    param_tab_layout->addWidget(label_lev, 0, 1, Qt::AlignCenter);
    param_tab_layout->addWidget(label_mean, 0, 2, Qt::AlignCenter);
    param_tab_layout->addWidget(label_uev, 0, 3, Qt::AlignCenter);

    param_tab_layout->addWidget(label_tau_inc, 1,0);
    param_tab_layout->addWidget(inc_lev, 1,1, Qt::AlignCenter);
    param_tab_layout->addWidget(inc_mean, 1,2, Qt::AlignCenter);
    param_tab_layout->addWidget(inc_uev, 1,3, Qt::AlignCenter);

    param_tab_layout->addWidget(label_percentage_predetect, 2,0);
    param_tab_layout->addWidget(percentage_predetection, 2,2, Qt::AlignCenter);

    param_tab_layout->addWidget(label_tau_symp, 3,0);
    param_tab_layout->addWidget(symp_lev, 3,1, Qt::AlignCenter);
    param_tab_layout->addWidget(symp_mean, 3,2, Qt::AlignCenter);
    param_tab_layout->addWidget(symp_uev, 3,3, Qt::AlignCenter);

    param_tab_layout->addWidget(label_asymp, 4, 0);
    param_tab_layout->addWidget(percentage_asymptomatic, 4, 2, Qt::AlignCenter);

    param_tab_layout->addWidget(label_tau_post, 5,0);
    param_tab_layout->addWidget(post_lev, 5,1, Qt::AlignCenter);
    param_tab_layout->addWidget(post_mean, 5,2, Qt::AlignCenter);
    param_tab_layout->addWidget(post_uev, 5,3, Qt::AlignCenter);

    param_tab_layout->addWidget(label_pcr_sens, 7, 0);
    param_tab_layout->addWidget(pcr_sens, 7, 2, Qt::AlignCenter);

    param_tab_layout->addWidget(label_pcr_spec, 8, 0);
    param_tab_layout->addWidget(pcr_spec, 8, 2, Qt::AlignCenter);

    param_tab_layout->addWidget(reset_PushButton, 9, 3);

    param_tab_layout->setHorizontalSpacing(20);
    param_tab_layout->setSizeConstraint(QLayout::SetFixedSize);

    QWidget *widget = new QWidget;
    widget->setLayout(param_tab_layout);
    return widget;
}

void MainWindow::initialize_test_date_checkboxes()
{
    // calc number of boxes needed
    int passed_days{0};
    passed_days = std::ceil(this->time_passed->value());
    int max_n_test = this->quarantine->value() + passed_days +1;

    // create boxes
    this->test_date_checkboxes.clear();
    QHBoxLayout *layout = new QHBoxLayout;
    for(int i=-passed_days; i < max_n_test - passed_days; ++i)
    {
        QCheckBox *box = new QCheckBox();
        if (i < 0 )
        {
            box->setEnabled(false);
        }
        else
        {
            this->test_date_checkboxes.push_back(box);
            if (i < int(this->test_date_checkboxes_states.size())){
                box->setChecked(this->test_date_checkboxes_states.at(i));
            }
        }

        QVBoxLayout *vbox = new QVBoxLayout;
        QLabel *label = new QLabel(QString::number(i));
        vbox->addWidget(box);
        vbox->addWidget(label);
        vbox->setAlignment(Qt::AlignTop);

        layout->addItem(vbox);
    }
    this->test_days_box->setLayout(layout);
    this->test_date_checkboxes_states.clear();
}

void MainWindow::update_test_date_checkboxes()
{
    //save states
    for (auto box : this->test_date_checkboxes)
    {
        if ( box->isEnabled() )
            this->test_date_checkboxes_states.push_back(box->isChecked());
    }
    //delete
    QLayout *hb = this->test_days_box->layout();
    while(!hb->isEmpty())
    {
        QLayout *vb = hb->takeAt(0)->layout();
        while(!vb->isEmpty())
        {
            QWidget *w = vb->takeAt(0)->widget();
            delete w;
        }
        delete vb;
    }
    delete hb;
    this->test_days_box->update();

    // re-initialize
    initialize_test_date_checkboxes();
}

void MainWindow::run_PushButton_clicked()
{
    Simulation *simulation = new Simulation(this);
    simulation->run();
}

void MainWindow::reset_PushButton_clicked()
{
    this->inc_lev->setValue(this->default_values["inc_lev"]);
    this->inc_mean->setValue(this->default_values["inc_mean"]);
    this->inc_uev->setValue(this->default_values["inc_uev"]);

    this->percentage_predetection->setValue(this->default_values["percentage_predetection"]);

    this->symp_lev->setValue(this->default_values["symp_lev"]);
    this->symp_mean->setValue(this->default_values["symp_mean"]);
    this->symp_uev->setValue(this->default_values["symp_uev"]);

    this->percentage_asymptomatic->setValue(this->default_values["percentage_asymptomatic"]);

    this->post_lev->setValue(this->default_values["post_lev"]);
    this->post_mean->setValue(this->default_values["post_mean"]);
    this->post_uev->setValue(this->default_values["post_uev"]);

    this->pcr_sens->setValue(this->default_values["pcr_sens"]);
    this->pcr_spec->setValue(this->default_values["pcr_spec"]);
}

void MainWindow::time_passed_valueChanged()
{
    update_test_date_checkboxes();
}

void MainWindow::quarantine_valueChanged()
{
    update_test_date_checkboxes();
}

void MainWindow::mode_ComboBox_currentIndexChanged(int)
{
    QLayout *input_tab_layout = this->tab->widget(0)->layout()->itemAt(0)->layout();
    QGridLayout *gridLayout = qobject_cast<QGridLayout*>(input_tab_layout);
    QWidget *widget = gridLayout->itemAtPosition(1,0)->widget();
    QLabel* label_time_passed = qobject_cast<QLabel*>(widget);

    label_time_passed->setText((std::string{"Time passed since "} +
                                this->mode_ComboBox->currentText().toStdString() +
                                std::string{" [days]:"}).c_str() );
}
