#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "simulation.h"

#include <numeric>
#include <QString>
#include <QLabel>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    set_checkboxes();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_checkBox_unknown_toggled()
{
    if ( ui->checkBox_unknown->isChecked() ){
        ui->input_timepassed->setEnabled(false);
        ui->checkBox_known->setChecked(false);
        ui->label_15->setEnabled(false);
    }
    else {
        ui->input_timepassed->setEnabled(true);
        ui->checkBox_known->setChecked(true);
        ui->label_15->setEnabled(true);
    }
    update_checkboxes();
}

void MainWindow::on_checkBox_known_toggled()
{
    if (ui->checkBox_known->isChecked()  ){
        ui->input_timepassed->setEnabled(true);
        ui->checkBox_unknown->setChecked(false);
        ui->label_15->setEnabled(true);
    }
    else {
        ui->input_timepassed->setEnabled(false);
        ui->checkBox_unknown->setChecked(true);
        ui->label_15->setEnabled(false);
    }
    update_checkboxes();
}

void MainWindow::set_checkboxes(){
    // calc number of boxes needed
    int passed_days{0};
    if( ui->checkBox_known->isChecked()){
        passed_days = std::ceil(ui->input_timepassed->value());
    }
    int max_n_test = ui->input_quarantaine->value() +passed_days +1;

    // resixe if necesarry (cb w14xh14)
    bool rescale = false;
    int rescale_size{};
    if (max_n_test*14 > ui->testDates->width()){
        rescale=true;
        rescale_size = std::floor(ui->testDates->width() / max_n_test);
    }

    // create boxes
    test_boxes.clear();
    QHBoxLayout *layout = new QHBoxLayout;
    for(int i=-passed_days; i < max_n_test-passed_days; ++i){
        QVBoxLayout *vbox = new QVBoxLayout;

        QCheckBox *box = new QCheckBox();
        if (rescale){
            box->setStyleSheet("QCheckBox::indicator { width:" +
                               QString::number(rescale_size) + "px; height:"+
                               QString::number(rescale_size) +
                               "px;}");
        }
        if (i < 0 ) { box->setEnabled(false);}
        else
        {
            test_boxes.push_back(box);
            if (i < int(test_boxes_states.size())){
                box->setChecked(test_boxes_states.at(i));
            }
        }

        QLabel *label = new QLabel(QString::number(i));
        vbox->addWidget(box);
        vbox->addWidget(label);

        layout->addItem(vbox);
    }
    ui->testDates->setLayout(layout);
    test_boxes_states.clear();
}
void MainWindow::update_checkboxes(){
    //save states
    for (auto box : test_boxes){
        if ( box->isEnabled() )
            test_boxes_states.push_back(box->isChecked());
    }
    //delete
    QLayout *hb = ui->testDates->layout();
    while(!hb->isEmpty()) {
        QLayout *vb = hb->takeAt(0)->layout();
        while(!vb->isEmpty()) {
            QWidget *w = vb->takeAt(0)->widget();
            delete w;
        }
        delete vb;
    }
    delete hb;
    ui->testDates->update();

    set_checkboxes();
}

void MainWindow::on_pushButton_run_clicked()
{
    Simulation *simulation = new Simulation(this);
    simulation->run();
}

void MainWindow::on_input_timepassed_valueChanged(double)
{
    update_checkboxes();
}

void MainWindow::on_input_quarantaine_valueChanged(double)
{
    update_checkboxes();
}
