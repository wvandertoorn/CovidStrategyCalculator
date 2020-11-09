#include "simulation.h"
#include "mainwindow.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <QVBoxLayout>
#include <QHeaderView>

#include <QtCharts/QAreaSeries>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QLegendMarker>
#include <QtCharts/QValueAxis>

std::vector<int> Simulation::sub_compartments = {5,1,13,1,1}; //predetect, presympt, symp, post, risk-node (sink is dropped)
int Simulation::nr_compartments = std::accumulate(sub_compartments.begin(), sub_compartments.end(), 0);
float Simulation::t_inf = 100.;

Simulation::Simulation(MainWindow *parent) : QObject(parent)
{
    this->m_parent = parent;
    collect_data(this->m_parent);

    initial_states.setZero(nr_compartments);

    switch (this->mode) {
        case 0: initial_states(0) = 1.0; // mode exposure
                break;
        case 1:
        {
                int first_symptomatic_compartment = sub_compartments[0] + sub_compartments[1];
                this->initial_states(first_symptomatic_compartment) = 1.0; // mode symptom onset
                initial_states(nr_compartments - 1) = 1.0;
                break;
        }
        case 2: break;
    }
}

Simulation::Simulation(MainWindow *parent, Eigen::MatrixXf init_states, float prerisk) : QObject(parent)
{
    this->m_parent = parent;
    collect_data(this->m_parent);

    init_states.resize(nr_compartments, 1);
    init_states(nr_compartments - 1, 0) = init_states(Eigen::seq(sub_compartments[0],
                                                                 sub_compartments[0]
                                                                 + sub_compartments[1]
                                                                 + sub_compartments[2] -1), 0)
                                                                .sum();

    this->initial_states = init_states;

    this->pre_test_infect_prob = prerisk;
}

std::tuple<std::vector<Eigen::MatrixXf>,
           std::vector<Eigen::MatrixXf>> Simulation::run_prevalence(std::vector<float> week_incidences)
{
    time_passed = 0;
    this->fraction_asymtomatic = 1.;
    this->t_test = {};
    //TODO document the fact that the prevalence estimator does not assume symptomatic screening in the first simulation.

    Eigen::MatrixXf states_today_mean, states_today_lev, states_today_uev,
                    probs_today_mean, probs_today_lev, probs_today_uev;

    states_today_mean.setZero(1, nr_compartments);
    states_today_lev.setZero(1, nr_compartments);
    states_today_uev.setZero(1, nr_compartments);

    float total_residence_mean = std::accumulate(residence_times_mean.begin(), residence_times_mean.end(), 0);
    // float tot_lev = std::accumulate(residence_times_lev.begin(), residence_times_lev.end(), 0);
    // float tot_uev = std::accumulate(residence_times_uev.begin(), residence_times_uev.end(), 0);

    for (int week = 0; week < 5; ++week)
    {
        float day_incidence = week_incidences[week] / 7;
        this->quarantine = (week + 1) * 7 - 1;

        // set intial states according to incidence report
        int counter = 0;
        for (int i = 0; i < 4; ++i)
        {
            for (int j=0; j < sub_compartments[i]; ++j)
            {
                    float probability = residence_times_mean[i] / total_residence_mean * day_incidence / sub_compartments[i];
                    initial_states(counter) = probability;
                    counter++;
            }
        }

        this->run();
        states_today_mean.array() += X_mean(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        states_today_lev.array() += X_lev(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        states_today_uev.array() += X_uev(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        //TODO think about how to handle risk node in prevalence estimator

        initial_states.array() = 0.0; //needed to set initial state of sink compartment back to zero
    }

    probs_today_mean = assemble_phases(states_today_mean, sub_compartments);
    probs_today_lev = assemble_phases(states_today_lev, sub_compartments);
    probs_today_uev = assemble_phases(states_today_uev, sub_compartments);

    std::vector<Eigen::MatrixXf> states{states_today_mean, states_today_lev, states_today_uev};
    std::vector<Eigen::MatrixXf> probs{probs_today_mean, probs_today_lev, probs_today_uev};

    return std::make_tuple(states, probs);
}

void Simulation::collect_data(MainWindow *parent)
{
    // strategy
    mode = parent->mode_ComboBox->currentIndex();
    mode_str = parent->mode_ComboBox->currentText().toStdString();
    time_passed = parent->time_passed->value();
    quarantine = parent->quarantine->value();
    use_symptomatic_screening = parent->use_symptomatic_screening->isChecked();

    // parameters
    if (use_symptomatic_screening){
        fraction_asymtomatic = parent->percentage_asymptomatic->value() /100;
    } else {fraction_asymtomatic = 1.;}

    // mean prediction
    residence_times_mean.clear();
    residence_times_mean.push_back(parent->percentage_predetection->value()/100 * parent->inc_mean->value());
    residence_times_mean.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_mean->value());
    residence_times_mean.push_back(parent->symp_mean->value());
    residence_times_mean.push_back(parent->post_mean->value());

    // lower extreme value prediction
    residence_times_lev.clear();
    residence_times_lev.push_back(parent->percentage_predetection->value()/100 * parent->inc_lev->value());
    residence_times_lev.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_lev->value());
    residence_times_lev.push_back(parent->symp_lev->value());
    residence_times_lev.push_back(parent->post_mean->value());

    // upper extreme value prediction
    residence_times_uev.clear();
    residence_times_uev.push_back(parent->percentage_predetection->value()/100 * parent->inc_uev->value());
    residence_times_uev.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_uev->value());
    residence_times_uev.push_back(parent->symp_uev->value());
    residence_times_uev.push_back(parent->post_mean->value());

    sensitivity = parent->pcr_sens->value() /100;
    specificity = parent->pcr_spec->value() /100;

    this->test_type = parent->test_type->currentText();

    switch (parent->test_type->currentIndex()) {
        case 1: sensitivity = sensitivity * (parent->rel_antigen_sens->value() /100);
                break; //antigen
    }

    t_test = collect_t_test(parent->test_date_checkboxes);
}

std::vector<int> Simulation::collect_t_test(std::vector<QCheckBox*> boxes)
{
    std::vector<int> v{};
    int counter = time_passed;
    for (auto box : boxes)
    {
        if ( box->isChecked() )
            v.push_back(counter);
        ++counter;
    }
    return v;
}

Eigen::VectorXf Simulation::calc_rates(std::vector<float> times,
                                       std::vector<int> comp
                                       )
{
    int n = std::accumulate(comp.begin(), comp.end(), 0);
    Eigen::VectorXf rates(n-1);
    int counter = 0;
    for (int i=0; i < 4; ++i)
    {
        for (int j=0; j < comp[i]; ++j)
        {
                rates[counter] = comp[i] /times[i];
                counter++;
        }
    }
    return rates;
}

Eigen::VectorXf Simulation::FOR_vector(std::vector<int> comp)
{
    int n = std::accumulate(comp.begin(), comp.end(), 0);
    Eigen::VectorXf rates;
    rates.setOnes(n);

    int counter = 0;
    for (int j=0; j < comp[0]; ++j)
    {
            rates[counter] = specificity;
            counter++;
    }
    for (int j=0; j < comp[1] + comp[2]; ++j)
    {
            rates[counter] = 1 - sensitivity;
            counter++;
    }

    return rates;
}

Eigen::MatrixXf Simulation::calc_S(int n)
{
    Eigen::MatrixXf S;
    S.setZero(n, n-1);

    for (int i=0; i < n-1; ++i)
    {
        S(i, i) = -1;
        S(i+1, i) = 1;
    }
    return S;
}

Eigen::MatrixXf Simulation::calc_A(Eigen::MatrixXf S,
                                   Eigen::VectorXf r,
                                   std::vector<int> comp)
{
    int n = S.rows();
    Eigen::MatrixXf temp;
    temp.setZero(n, n-1);

    for (int i = 0; i < n;  ++i)
    {
        for (int j=0; j < n-1; ++j)
        {
            temp(i,j) = r(j);
        }
    }

    Eigen::MatrixXf A;
    A = S.cwiseProduct(temp);

    Eigen::MatrixXf A_square(n-1, n-1);
    A_square = A(Eigen::seq(0, Eigen::last-1), Eigen::all); //drop last row, sink node not of interest

    Eigen::MatrixXf A_augmented;
    A_augmented.setZero(n,n);
    A_augmented(Eigen::seq(0, Eigen::last-1),Eigen::seq(0, Eigen::last-1)) = A_square;
    A_augmented(Eigen::last, Eigen::seq(comp[0], comp[0] + comp[1] - 1)).array() = 1.;
    A_augmented(Eigen::last, Eigen::seq(comp[0] + comp[1], comp[0] + comp[1] + comp[2] - 1)).array() = fraction_asymtomatic;

    return A_augmented;
}

Eigen::MatrixXf  Simulation::calc_X(int delay,
                                    int qrntn,
                                    Eigen::MatrixXf A,
                                    Eigen::VectorXf states
                                    )
{
    int t_end = delay + qrntn + 1;

    Eigen::VectorXf time(t_end);
    for (int i=0; i< t_end; ++i)
    {
        time(i) = float(i);
    }

    int n = states.size();
    Eigen::MatrixXf X;
    X.setZero(t_end, n);

    for (int i=0; i<t_end; ++i)
    {
        X.row(i) =(A * time[i]).exp()*states;
    }

    return X;
}

float Simulation::calc_risk_at_T(Eigen::MatrixXf A,
                                Eigen::VectorXf states,
                                float time_T)
{
    int n = states.size();
    Eigen::MatrixXf X;
    X.setZero(1, n);

    X.row(0) = (A * time_T).exp()*states;

    return X(0, Eigen::last);
}

Eigen::MatrixXf Simulation::assemble_phases(Eigen::MatrixXf X,
                                            std::vector<int> comp
                                            )
{
    int n_time = X.rows();
    Eigen::MatrixXf assembled(n_time, 5);
    int col_counter = 0;
    for (int i=0; i<5; ++i)
    {
        assembled(Eigen::all, i) = X(Eigen::all, Eigen::seq(col_counter,
                                                            col_counter + comp.at(i) -1))
                                    .rowwise()
                                    .sum();
        col_counter = col_counter + comp.at(i);
    }
    return assembled;
}

Eigen::MatrixXf Simulation::risk_node_to_relative_residual_risk(Eigen::MatrixXf risk, float risk_T){

    Eigen::MatrixXf risk_T_vector(risk.rows(), 1);
    risk_T_vector.fill(risk_T);

    risk.array() = risk_T_vector.array() - risk.array();

    risk.array() = (risk.array() / risk(0,0)).array() * 100.;
    return risk;
}


QtCharts::QChartView* Simulation::create_plot(Eigen::MatrixXf detectable,
                                              Eigen::MatrixXf mean,
                                              Eigen::MatrixXf uev,
                                              Eigen::MatrixXf lev,
                                              std::vector<int> time_range_for_plot
                                        )
{
    QtCharts::QChart *chart = new QtCharts::QChart();

    QtCharts::QLineSeries *detect_low = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *detect_high = new QtCharts::QLineSeries;

    for (int j=0; j < detectable.rows(); ++j)
    {
        detect_low->append(j-time_passed, detectable(j,0));
        detect_high->append(j-time_passed, detectable(j,1));
    }
    QtCharts::QAreaSeries *detect_area = new QtCharts::QAreaSeries(detect_low, detect_high);
    detect_area->setName("Assay sensitivity");

    QColor detect_color(140, 129, 152, 128);
    QPen detect_pen(detect_color);
    detect_area->setPen(detect_pen);
    QBrush detect_brush(detect_color);
    detect_area->setBrush(detect_brush);

    chart->addSeries(detect_area);

    QtCharts::QLineSeries *risk_low = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *risk_mean = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *risk_high = new QtCharts::QLineSeries;
    risk_mean->setName("Relative residual risk");

    for (int j = 0; j <  (int) time_range_for_plot.size(); ++j)
    {
        risk_low->append(time_range_for_plot[j], lev(j, 0));
        risk_mean->append(time_range_for_plot[j], mean(j, 0));
        risk_high->append(time_range_for_plot[j], uev(j, 0));
    }
    QtCharts::QAreaSeries *risk_area = new QtCharts::QAreaSeries(risk_low, risk_high);

    QColor risk_color(255, 128, 128, 128);
    QPen risk_pen(risk_color);
    risk_area->setPen(risk_pen);
    QBrush risk_brush(risk_color);
    risk_area->setBrush(risk_brush);
    QPen red_pen(Qt::red);
    risk_mean->setPen(red_pen);

    chart->addSeries(risk_area);
    chart->addSeries(risk_mean);
    chart->legend()->markers(risk_area)[0]->setVisible(false);

    //------------------

    chart->createDefaultAxes();
    switch (this->mode) {
        case 0: chart->axes(Qt::Horizontal).first()->setTitleText("Days since infection");
                break;
        case 1: chart->axes(Qt::Horizontal).first()->setTitleText("Days since symptom onset");
                break;
        case 2: chart->axes(Qt::Horizontal).first()->setTitleText("Days since exposure");
                break;
    }
    chart->axes(Qt::Vertical).first()->setTitleText("Percent");

    QtCharts::QValueAxis *axisX = qobject_cast<QtCharts::QValueAxis*>(chart->axes(Qt::Horizontal).first());
    axisX->setTickCount(time_passed + quarantine  + 1);

    QtCharts::QValueAxis *axisY = qobject_cast<QtCharts::QValueAxis*>(chart->axes(Qt::Vertical).first());
    axisY->setRange(0, 100);

    chart->legend()->show();
    QFont f("Helvetica", 10);
    chart->legend()->setFont(f);

    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    return chartView;
}

QTableWidget* Simulation::create_table(Eigen::MatrixXf detectable,
                                       Eigen::MatrixXf mean,
                                       Eigen::MatrixXf uev,
                                       Eigen::MatrixXf lev)
{
    int n_time = detectable.rows();
    int offset = time_passed;

    QStringList h_labels;
    for (int i = - offset; i < n_time - offset; ++i)
    {
        h_labels << QString::number(i);
    }

    QTableWidget *table = new QTableWidget(1, n_time);
    table->setHorizontalHeaderLabels(h_labels);
    table->setVerticalHeaderLabels((QStringList() << "Assay sensitivity [%]" << "Relative residual risk [%]"));

    for (int i = 0; i< n_time; ++i)
    {
        table->setItem(0, i, new QTableWidgetItem("("+ QString::number(detectable(i, 0), 'f', 2)
                                                  + ", " + QString::number(detectable(i, 1), 'f', 2) + ")"));
        // table->setItem(1, i, new QTableWidgetItem(QString::number(mean(i, 0), 'f', 2) +
        //                                           + "  (" + QString::number(lev(i, 0), 'f', 2)
        //                                           + ", " + QString::number(uev(i, 0), 'f', 2) + ")"));

        table->item(0, i)->setFlags(table->item(0, i)->flags() &  ~Qt::ItemIsEditable);
        // table->item(1, i)->setFlags(table->item(1, i)->flags() &  ~Qt::ItemIsEditable);
    }
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    return table;
}

void Simulation::create_result_log()
{
    QLabel *label = new QLabel(tr("Result log"));

    QTableWidget *table = new QTableWidget(1, 8);
    table->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    table->setHorizontalHeaderLabels((QStringList() << "mode"
                                                    << "sympt screening \n (perc asympt)"
                                                    << "time passed \n [days]"
                                                    << "quarantine/isolation \n [days]"
                                                    << "test \n [days]"
                                                    << "test type"
                                                    << "start risk \n [%]"
                                                    << "fold risk \n reduction"));

    write_row_result_log(table);
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    connect(
        table->horizontalHeader(),
        SIGNAL(sectionResized(int, int, int)),
        table,
        SLOT(resizeRowsToContents()));

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(label, 0);
    layout->addWidget(table, 90);

    this->m_parent->log->setLayout(layout);
}

void Simulation::write_row_result_log(QTableWidget *table)
{
    table->setItem(0, 0, new QTableWidgetItem(QString(this->mode_str.c_str())));

    QString boolText = this->use_symptomatic_screening ? "yes" : "no";

    if (this->use_symptomatic_screening){
        table->setItem(0, 1, new QTableWidgetItem(boolText
                                                 + " ("
                                                 + QString::number(this->fraction_asymtomatic * 100)
                                                 + "%)"));
    } else {table->setItem(0, 1, new QTableWidgetItem(boolText));}

    table->setItem(0, 2, new QTableWidgetItem(QString::number(this->time_passed)));
    table->setItem(0, 3, new QTableWidgetItem(QString::number(this->quarantine)));

    if (this->t_test.size())
    {
        QString days{};
        for (int day : this->t_test)
        {
            days += (QString::number(day - this->time_passed) + ", ");
        }
        days.chop(2);
        table->setItem(0, 4, new QTableWidgetItem(days));

        table->setItem(0, 5, new QTableWidgetItem(this->test_type));
    }
    else
    {
      table->setItem(0, 4, new QTableWidgetItem());
      table->setItem(0, 5, new QTableWidgetItem());
    }

    table->setItem(0, 6, new QTableWidgetItem(QString::number(pre_test_infect_prob *100, 'f', 2)));

    // table->setItem(0, 7, new QTableWidgetItem(QString::number(calculate_strategy_result(result_matrix_mean)* 100,
    //                                                          'f', 2)
    //                                          + "  ("
    //                                          + QString::number(calculate_strategy_result(result_matrix_uev)* 100,
    //                                                            'f', 2)
    //                                          + ", "
    //                                          + QString::number(calculate_strategy_result(result_matrix_lev)* 100,
    //                                                            'f', 2)
    //                                          + ")" ));
    table->setItem(0, 7, new QTableWidgetItem(QString::number(this->fold_RR_mean,
                                                             'f', 2)
                                             + "  ("
                                             + QString::number(this->fold_RR_uev,
                                                               'f', 2)
                                             + ", "
                                             + QString::number(this->fold_RR_lev,
                                                               'f', 2)
                                             + ")" ));

    for (int i = 0; i < 8; ++i)
    {
        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
}

void Simulation::update_result_log()
{
    QWidget *widget = this->m_parent->log->layout()->takeAt(1)->widget();
    QTableWidget* table = qobject_cast<QTableWidget*>(widget);
    table->insertRow(0);
    write_row_result_log(table);

    this->m_parent->log->layout()->addWidget(table);
}

float Simulation::calculate_strategy_result(Eigen::MatrixXf matrix)
{
    float risk = 1.;

    if (t_test.size())
    {
        for (int day : t_test)
        {
            risk = risk * (1. - sensitivity) * (matrix(day, 1) + fraction_asymtomatic * matrix(day, 2) ) +
                   risk * specificity * matrix(day, 0);
        }
        if (t_test.back() < matrix.rows()-1)
        {
            risk = risk * (matrix(Eigen::last, Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(Eigen::last, 2))
                    / (matrix(t_test.back(), Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(t_test.back(), 2));
        }
    }
    else risk = matrix(Eigen::last, Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(Eigen::last, 2);
    return risk;
}

Eigen::MatrixXf Simulation::calculate_assay_sensitivity()
{
    int n_time = assay_detectibility_mean_case.rows();
    Eigen::MatrixXf detectibility(n_time, 2);

    for (int j=0; j<n_time; ++j)
    {
        float m = assay_detectibility_mean_case(j, 0) * (1 - specificity) + assay_detectibility_mean_case(j, Eigen::seq(1, 3)).sum() * sensitivity;
        float lev = assay_detectibility_worst_case(j, 0) * (1 - specificity) + assay_detectibility_worst_case(j, Eigen::seq(1, 3)).sum() * sensitivity;
        float uev = assay_detectibility_best_case(j, 0) * (1 - specificity) + assay_detectibility_best_case(j, Eigen::seq(1, 3)).sum() * sensitivity;

        std::vector<float> v{m, lev, uev};

        int maxElementIndex = std::max_element(v.begin(), v.end()) - v.begin();
        int minElementIndex = std::min_element(v.begin(), v.end()) - v.begin();

        detectibility(j, 0) =  v[minElementIndex];
        detectibility(j, 1) =  v[maxElementIndex];
    }

    float infected_population = initial_states(Eigen::seq(0, Eigen::last - 1)).sum();
    return detectibility.array() / infected_population * 100.;
}

void Simulation::output_results()
{
    Eigen::MatrixXf mean = risk_node_to_relative_residual_risk( X_mean(Eigen::all, Eigen::last), risk_T_preprocedure_mean);
    Eigen::MatrixXf uev = risk_node_to_relative_residual_risk( X_uev(Eigen::all, Eigen::last), risk_T_preprocedure_uev);
    Eigen::MatrixXf lev = risk_node_to_relative_residual_risk( X_lev(Eigen::all, Eigen::last), risk_T_preprocedure_lev);

    Eigen::MatrixXf detectable = calculate_assay_sensitivity();

    QtCharts::QChartView* plot = create_plot(detectable, mean, uev, lev, time_for_plot);

    QTableWidget* table = create_table(detectable, mean, uev, lev);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setSpacing( 0 );

    layout->addWidget(plot, 85);
    layout->addWidget(table, 8);

    // if previous layout exists, delete
    if (this->m_parent->chart->layout())
    {
        QLayout *hb = this->m_parent->chart->layout();
        while(!hb->isEmpty())
        {
            QWidget *w = hb->takeAt(0)->widget();
            delete w;
        }
        delete hb;
    }

    this->m_parent->chart->setLayout(layout);
    this->m_parent->chart->update();

    // initialize or update
    if (!this->m_parent->log->layout())
    {
        create_result_log();
    }
    else
    {
        update_result_log();
    }

}

void Simulation::run()
{
    Eigen::VectorXf rates;
    Eigen::MatrixXf S;
    Eigen::MatrixXf A;
    Eigen::MatrixXf X;
    float pre_procedure, post_procedure;

    S = calc_S(nr_compartments);

    std::vector<float> worst_case_extreme{residence_times_uev[0],
                                          residence_times_lev[1],
                                          residence_times_lev[2],
                                          residence_times_lev[3]};
    std::vector<float> best_case_extreme{residence_times_lev[0],
                                         residence_times_uev[1],
                                         residence_times_uev[2],
                                         residence_times_uev[3]};

    // mean case
    rates = calc_rates(residence_times_mean, sub_compartments);
    A = calc_A(S, rates, sub_compartments);
    if (t_test.size()){
        auto results = calculate_strategy_with_test(A);
        this->X_mean = std::get<0>(results);
        this->result_matrix_mean = assemble_phases(X_mean, sub_compartments);
        this->time_for_plot = std::get<1>(results);
        this->risk_T_postprocedure_mean = std::get<2>(results);
        this->risk_T_preprocedure_mean = calc_risk_at_T(A, initial_states,  t_inf);

        pre_procedure =  risk_T_preprocedure_mean - initial_states(Eigen::last);
        post_procedure = risk_T_postprocedure_mean - X_mean(Eigen::last, Eigen::last);
    }
    else {
        this->X_mean = calc_X(time_passed, quarantine, A, initial_states);
        this->result_matrix_mean = assemble_phases(X_mean, sub_compartments);
        this->risk_T_postprocedure_mean = calc_risk_at_T(A, initial_states,  t_inf);
        this->risk_T_preprocedure_mean = risk_T_postprocedure_mean;

        pre_procedure = risk_T_preprocedure_mean - X_mean(0, Eigen::last);
        post_procedure = risk_T_postprocedure_mean - X_mean(Eigen::last, Eigen::last);

        std::vector<int> v(time_passed + quarantine + 1);
        std::iota(v.begin(), v.end(), -time_passed);
        this->time_for_plot = v;
    }
    this->fold_RR_mean = pre_procedure / post_procedure;

    X = calc_X(time_passed, quarantine, A, initial_states);
    this->assay_detectibility_mean_case = assemble_phases(X, sub_compartments);

    // worst case in terms of detectibility (long predetect, short rest )
    rates = calc_rates(worst_case_extreme, sub_compartments);
    A = calc_A(S, rates, sub_compartments);
    if (t_test.size()){
        auto results = calculate_strategy_with_test(A);
        this->X_lev = std::get<0>(results);
        this->result_matrix_lev = assemble_phases(X_lev, sub_compartments);
        this->risk_T_postprocedure_lev = std::get<2>(results);
        this->risk_T_preprocedure_lev = calc_risk_at_T(A, initial_states,  t_inf);

        pre_procedure = risk_T_preprocedure_lev - initial_states(Eigen::last);
        post_procedure = risk_T_postprocedure_lev - X_lev(Eigen::last, Eigen::last);

    }
    else {
        this->X_lev = calc_X(time_passed, quarantine, A, initial_states);
        this->result_matrix_lev = assemble_phases(X_lev, sub_compartments);
        this->risk_T_postprocedure_lev = calc_risk_at_T(A, initial_states,  t_inf);
        this->risk_T_preprocedure_lev = risk_T_postprocedure_lev;

        pre_procedure = risk_T_preprocedure_lev - X_lev(0, Eigen::last);
        post_procedure = risk_T_postprocedure_lev - X_lev(Eigen::last, Eigen::last);
    }
    this->fold_RR_lev = pre_procedure / post_procedure;

    X = calc_X(time_passed, quarantine, A, initial_states);
    this->assay_detectibility_worst_case = assemble_phases(X, sub_compartments);

    // best case in terms of detectibility (short predetect, long rest )
    rates = calc_rates(best_case_extreme, sub_compartments);
    A = calc_A(S, rates, sub_compartments);
    if (t_test.size()){
        auto results = calculate_strategy_with_test(A);
        this->X_uev = std::get<0>(results);
        this->result_matrix_uev = assemble_phases(X_uev, sub_compartments);
        this->risk_T_postprocedure_uev = std::get<2>(results);
        this->risk_T_preprocedure_uev = calc_risk_at_T(A, initial_states,  t_inf);

        pre_procedure =  risk_T_preprocedure_uev - initial_states(Eigen::last);
        post_procedure = risk_T_postprocedure_uev - X_uev(Eigen::last, Eigen::last);
    }
    else {
        this->X_uev = calc_X(time_passed, quarantine, A, initial_states);
        this->result_matrix_uev = assemble_phases(X_uev, sub_compartments);
        this->risk_T_postprocedure_uev = calc_risk_at_T(A, initial_states,  t_inf);
        this->risk_T_preprocedure_uev = risk_T_postprocedure_uev;

        pre_procedure = risk_T_preprocedure_uev - X_uev(0, Eigen::last);
        post_procedure = risk_T_postprocedure_uev - X_uev(Eigen::last, Eigen::last);
    }
    this->fold_RR_uev = pre_procedure / post_procedure;

    X = calc_X(time_passed, quarantine, A, initial_states);
    this->assay_detectibility_best_case = assemble_phases(X, sub_compartments);

    calculate_assay_sensitivity();
}

std::tuple<Eigen::MatrixXf, std::vector<int>, float> Simulation::calculate_strategy_with_test(Eigen::MatrixXf A){

    Eigen::MatrixXf X0(nr_compartments, 1);
    X0.array() = this->initial_states.array();

    Eigen::MatrixXf X_total(this->quarantine + this->time_passed + 1 + this->t_test.size(), nr_compartments);

    X0.resize(1, nr_compartments);
    X_total(0, Eigen::all).array() = X0.array();
    X0.resize(nr_compartments, 1);

    Eigen::MatrixXf FOR(nr_compartments, 1);
    FOR.array() = FOR_vector(sub_compartments).array();
    FOR.resize(1, nr_compartments);

    std::vector<int> time_plot{ 0 - this->time_passed};

    int days_till_next_test{0}, days_till_end{0};
    int row_counter = 1;
    int t_counter = time_plot.back();

    // test_day is counted relative from exposure date, i.e. test at day 3 with delay of 3 days: test_day==6
    for (int test_day : this->t_test){
        test_day = test_day - this->time_passed;
        days_till_next_test = test_day - t_counter;

        std::vector<int> v(days_till_next_test);
        std::iota(v.begin(), v.end(), t_counter + 1);
        time_plot.insert( time_plot.end(), v.begin(), v.end() );
        time_plot.push_back(t_counter + days_till_next_test);

        Eigen::MatrixXf X = calc_X(days_till_next_test, 0, A, X0);

        Eigen::MatrixXf final_state(1, nr_compartments);
        final_state.array() = X(Eigen::last, Eigen::all).array() * FOR.array();

        final_state(0, nr_compartments - 1) += (X(Eigen::last, Eigen::seq(sub_compartments[0], Eigen::last)).array()
                                                - final_state(0, Eigen::seq(sub_compartments[0], Eigen::last)).array())
                                               .sum();

        X_total(Eigen::seq(row_counter, row_counter + days_till_next_test -1), Eigen::all).array() = X(Eigen::seq(1, Eigen::last), Eigen::all).array();
        X_total(row_counter + days_till_next_test , Eigen::all).array() = final_state.array();
        final_state.resize(nr_compartments, 1);
        X0.array() = final_state.array();
        row_counter = time_plot.size();
        t_counter = time_plot.back();
    }
    if (this->quarantine > t_counter){
        days_till_end = this->quarantine - t_counter;
        Eigen::MatrixXf X = calc_X(days_till_end, 0, A, X0);
        X_total(Eigen::seq(row_counter, row_counter + days_till_end -1), Eigen::all).array() = X(Eigen::seq(1, Eigen::last), Eigen::all).array();

        std::vector<int> v(days_till_end);
        std::iota(v.begin(), v.end(), t_counter + 1);
        time_plot.insert( time_plot.end(), v.begin(), v.end() );
    }

    Eigen::MatrixXf tmp = X_total(Eigen::last, Eigen::all);
    tmp.resize(nr_compartments, 1);
    X0.array() = tmp.array();

    float risk_T = calc_risk_at_T( A, X0, this->t_inf - (float) this->quarantine - (float) this->time_passed);

    return std::make_tuple(X_total, time_plot, risk_T);
}
