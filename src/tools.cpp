#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
    const vector<VectorXd> &ground_truth) 
{
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
        || estimations.size() == 0){
            cout << "Invalid estimation or ground_truth data" << endl;
            return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i)
    {

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

void Tools::EstimateMeasurementNoise(const vector<MeasurementPackage>& measurements, const vector<GroundTruthPackage>& ground_truth)
{
    const size_t N = measurements.size();
    if (ground_truth.size() != N || ground_truth.size() == 0)
    {
        cout << "Size error. Cannot proceed!" << endl;
        return;
    }

    // Calculate variances via Welford's method
    double m_x_laser = 0.0;
    double s_x_laser = 0.0;
    double m_y_laser = 0.0;
    double s_y_laser = 0.0;

    double m_ro_radar      = 0.0;
    double s_rho_radar     = 0.0;
    double m_theta_radar   = 0.0;
    double s_theta_radar   = 0.0;
    double m_ro_dot_radar  = 0.0;
    double s_ro_dot_radar = 0.0;

    long N_laser           = 0;
    long N_radar           = 0;

    for (size_t k = 0; k < N; ++k) {
        const double x_gt  =  ground_truth[k].gt_values_(0);
        const double y_gt  =  ground_truth[k].gt_values_(1);
        const double vx_gt =  ground_truth[k].gt_values_(2);
        const double vy_gt =  ground_truth[k].gt_values_(3);

        if (measurements[k].sensor_type_ == MeasurementPackage::LASER)
        {
            const double x_meas       = measurements[k].raw_measurements_(0);
            const double y_meas       = measurements[k].raw_measurements_(1); 

            const double innovation_x = x_gt - x_meas;
            const double innovation_y = y_gt - y_meas;

            const double oldm_x_laser = m_x_laser;
            const double oldm_y_laser = m_y_laser;

            m_x_laser                += (innovation_x-m_x_laser)/(k+1);
            s_x_laser                += (innovation_x-m_x_laser)*(innovation_x-oldm_x_laser);

            m_y_laser                += (innovation_y-m_y_laser)/(k+1);
            s_y_laser                += (innovation_y-m_y_laser)*(innovation_y-oldm_y_laser);

            ++N_laser;
        }
        else if (measurements[k].sensor_type_ == MeasurementPackage::RADAR) {
            // output the estimation in the cartesian coordinates
            const double ro_gt     = sqrt(x_gt*x_gt+y_gt*y_gt);
            const double theta_gt  = atan2(y_gt,x_gt);

            if (fabs(ro_gt)>0.0001)
            {
                const double ro_dot_gt = (x_gt*vx_gt+y_gt*vy_gt)/ro_gt;

                const double ro_meas      = measurements[k].raw_measurements_(0);
                const double theta_meas   = measurements[k].raw_measurements_(1); 
                const double ro_dot_meas  = measurements[k].raw_measurements_(2);

                const double innovation_ro     = ro_gt - ro_meas;
                const double innovation_theta  = theta_gt - theta_meas;
                const double innovation_ro_dot = ro_dot_gt - ro_dot_meas;

                const double oldm_ro_radar = m_ro_radar;
                const double oldm_theta_radar = m_theta_radar;
                const double oldm_ro_dot_radar = m_ro_dot_radar;

                m_ro_radar     += (innovation_ro-m_ro_radar)/(k+1);
                s_rho_radar    += (innovation_ro-m_ro_radar)*(innovation_ro-oldm_ro_radar);
                m_theta_radar  += (innovation_theta-m_theta_radar)/(k+1);
                s_theta_radar  += (innovation_theta-m_theta_radar)*(innovation_theta-oldm_theta_radar);
                m_ro_dot_radar += (innovation_ro_dot-m_ro_dot_radar)/(k+1);
                s_ro_dot_radar += (innovation_ro_dot-m_ro_dot_radar)*(innovation_ro_dot-oldm_ro_dot_radar);

                ++N_radar;
            }
        }
        else
        {
            continue;
        }
    }

    if (N_laser > 1)
    {
        const double var_x_laser = s_x_laser/(N_laser-1);
        const double var_y_laser = s_y_laser/(N_laser-1);
        cout << "VarX Laser: " << var_x_laser << "VarY Laser: " << var_y_laser << endl;
    }

    if (N_radar > 1)
    {
        const double var_rho_radar     = s_rho_radar/(N_radar-1);
        const double var_theta_radar   = s_theta_radar/(N_radar-1);
        const double var_ro_dot_radar = s_ro_dot_radar/(N_radar-1);
        cout << "Var_Ro Radar: " << var_rho_radar << "Var_Theta Radar: " << var_theta_radar << "Var_Ro_Dot Radar: " << var_ro_dot_radar << endl;
    }
}

void Tools::EstimateProcessNoise(const vector<GroundTruthPackage>& ground_truth)
{
    const size_t N  = ground_truth.size();
    double m_a      = 0.0;
    double s_a      = 0.0;

    double m_psidd  = 0.0;
    double s_psidd  = 0.0;

    long long num_a      = 0;
    long long num_psi_dd = 0;

    for (size_t k = 1; k < N; ++k) {
        const double vx_gt     =  ground_truth[k].gt_values_(2);
        const double vy_gt     =  ground_truth[k].gt_values_(3);
        const double v_gt      =  sqrt(vx_gt*vx_gt+vy_gt*vy_gt);
        const long long t1     =  ground_truth[k].timestamp_;

        const double vx_old_gt =  ground_truth[k-1].gt_values_(2);
        const double vy_old_gt =  ground_truth[k-1].gt_values_(3);
        const double v_old_gt  =  sqrt(vx_old_gt*vx_old_gt+vy_old_gt*vy_old_gt);
        const long long t0     =  ground_truth[k-1].timestamp_;

        const double dt        = (t1-t0)/ 1000000.0;

        if (dt > 0.00001)
        {
            const double a = (v_gt - v_old_gt)/dt;

            const double m_a_old = m_a;
            m_a += (a-m_a)/(num_a+1);
            s_a += (a-m_a)*(a-m_a_old);

            ++num_a;

            if (k > 1)
            {
                const double psi_next = atan2(vy_gt, vx_gt);
                const double psi_curr = atan2(vy_old_gt, vx_old_gt);
                const double psi_last = atan2(ground_truth[k-2].gt_values_(3), ground_truth[k-2].gt_values_(2));

                double psi_diff = psi_next - psi_curr;
                //angle normalization
                while (psi_diff> M_PI) psi_diff-=2.*M_PI;
                while (psi_diff<-M_PI) psi_diff+=2.*M_PI;

                double last_psi_diff = psi_curr - psi_last;
                //angle normalization
                while (last_psi_diff> M_PI) last_psi_diff-=2.*M_PI;
                while (last_psi_diff<-M_PI) last_psi_diff+=2.*M_PI;

                const double psidd    = (psi_diff - last_psi_diff)/(dt*dt);

                const double m_psidd_old = m_psidd;
                m_psidd += (psidd-m_psidd)/(num_psi_dd+1);
                s_psidd += (psidd-m_psidd)*(psidd-m_psidd_old);

                ++num_psi_dd;
            }
        }    
    }

    if(num_a > 1)
    {
        const double var_a = s_a/(num_a-1);
        cout << "Q for A: " << var_a << endl;
    }

    if (num_psi_dd > 1)
    {
        const double var_psidd = s_psidd/(num_psi_dd-1);
        cout << "Q for PsiDD: " << var_psidd << endl;
    }

}