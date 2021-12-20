#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H


//Assignment 3 code 
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <init_state.h>
#include <find_min_vertices.h>
#include <fixed_point_constraints.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <T_linear_tetrahedron.h>
#include <V_linear_tetrahedron.h>
#include <V_spring_particle_particle.h>
#include <dV_linear_tetrahedron_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <mass_matrix_mesh.h>
#include <assemble_forces.h>
#include <assemble_stiffness.h>
#include <linearly_implicit_euler.h>
#include <implicit_euler.h>
#include <build_skinning_matrix.h>

//Variable for geometry
Eigen::MatrixXd V; //vertices of simulation mesh 
Eigen::MatrixXi T; //faces of simulation mesh
Eigen::MatrixXi F; //faces of simulation mesh

//variables for skinning
Eigen::MatrixXd V_skin;
Eigen::MatrixXi F_skin;
Eigen::SparseMatrixd N; 

//material parameters
double density = 0.1;
double YM = 6e5; //young's modulus
double mu = 0.4; //poissons ratio
double D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
double C = 0.5*YM/(2.0*(1.0+mu));

//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::SparseMatrixd P;
Eigen::VectorXd x0; 

//mass matrix
Eigen::SparseMatrixd M;
Eigen::VectorXd v0;

//scratch memory for assembly
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::SparseMatrixd tmp_stiffness;

using std::endl; using std::string;
using std::map; using std::copy;
std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points;

bool skinning_on = true;
bool fully_implicit = false;
bool bunny = true;
bool simulate_vol = false;
const char *energies[5] = { "smith_14", "smith_13",
                          "bower", "wang","ogden" };
std::vector<unsigned int> constant_pickeds;
float oom_ke = 0;
float oom_pe = 0;
int initial = 0;
int cur_energy = 0;
std::string energy_type=energies[cur_energy];
map<string, float> range_ke_map = {{"bower_arma",0.0001},
                                   {"smith_14_arma",0.0001},
                                   {"smith_13_arma",0.0001},
                                   {"ogden_arma",0.001},
                                   {"ogden_bunny",1e6},
                                   {"smith_13_bunny",1e6},
                                   {"smith_14_bunny",1e7},
                                   {"bower_bunny",1e4},
};
map<string, float> range_pe_map = {{"bower_arma",1},
                                   {"smith_14_arma",1000},
                                   {"smith_13_arma",1},
                                   {"ogden_arma",1},
                                   {"ogden_bunny",1e7},
                                   {"smith_13_bunny",1e7},
                                   {"smith_14_bunny",1e12},
                                   {"bower_bunny",1e4},
                         };



//selection spring
double k_selected = 1e5;

//force multiplier
double fm = 1.0;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {

    double V_ele, T_ele, KE,PE;

    spring_points.clear();

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
    
    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() * fm + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));

    }

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> qdot_1)->double {
        double E = 0;
        Eigen::VectorXd newq = P.transpose()*(q+dt*qdot_1)+x0;

        for(unsigned int ei=0; ei<T.rows(); ++ei) {
            
            V_linear_tetrahedron(V_ele,newq , V, T.row(ei), v0(ei), C, D, energy_type);
            E += V_ele;
        }

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {   
            V_spring_particle_particle(V_ele, spring_points[pickedi].first, newq.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            E += V_ele;
        }

        E += 0.5*(qdot_1 - qdot).transpose()*M*(qdot_1 - qdot);

        return E;
    };

    auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) { 
        
            assemble_forces(f, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D, energy_type);

            for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
                dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose()*q2+x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
            }

            f = P*f; 
        };

    //assemble stiffness matrix,
    auto stiffness = [&](Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) {
        assemble_stiffness(K, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D, energy_type);
        K = P*K*P.transpose();
    };

        if(fully_implicit)
            implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_qdot, tmp_force, tmp_stiffness);
        else
            linearly_implicit_euler(q, qdot, dt, M, force, stiffness, tmp_force, tmp_stiffness);
        

        
    KE = 0;
    PE = 0;

    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D, energy_type);
        PE += V_ele;
    }
    if(initial < 10 or pow(10,2*oom_ke)< KE or pow(10,2*oom_pe)< PE){
        oom_ke = floor(log10(abs(KE)));
        oom_pe = floor(log10(abs(PE)));
        initial +=1;
    }
    string bun = "";
    if(bunny){
       bun = "bunny";
    }else{
        bun="arma";
    }
    KE = KE / ((range_ke_map[energy_type +"_"+bun])*fm);
    PE = PE / (range_pe_map[energy_type +"_"+bun]*fm);
    Visualize::add_energy(t, KE, PE);
//     std::cout << t << "\t" << KE << "\t" << PE << "\t"<< oom_pe<<std::endl;
}
double calc_volume(double &energy, Eigen::Ref<const Eigen::VectorXd> q,
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element){

    Eigen::Vector3d a = q.segment<3>(3 * element(1)) - q.segment<3>(3 * element(0));
    Eigen::Vector3d b = q.segment<3>(3 * element(2)) - q.segment<3>(3 * element(0));
    Eigen::Vector3d c = q.segment<3>(3 * element(3)) - q.segment<3>(3 * element(0));
    energy = (a.cross(b)).dot(c)/6.0;
}
inline void simulate_volume(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {

    double V_ele, T_ele, KE,PE,VG, VG_ele;

    spring_points.clear();


    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
    Eigen::Vector3d constant_f;
    constant_f << 0.002813, 0.002813, 0;
    if(!Visualize::picked_vertices().empty()){
        constant_pickeds.push_back(Visualize::picked_vertices()[0]);

    }
    for(unsigned int pickedi = 0; pickedi < constant_pickeds.size(); pickedi++) {
        spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*constant_pickeds[pickedi]) + constant_f * fm + Eigen::Vector3d::Constant(1e-6),3*constant_pickeds[pickedi]));

    }

    auto energy = [&](Eigen::Ref<const Eigen::VectorXd> qdot_1)->double {
        double E = 0;
        Eigen::VectorXd newq = P.transpose()*(q+dt*qdot_1)+x0;

        for(unsigned int ei=0; ei<T.rows(); ++ei) {

            V_linear_tetrahedron(V_ele,newq , V, T.row(ei), v0(ei), C, D, energy_type);
            E += V_ele;
        }

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
            V_spring_particle_particle(V_ele, spring_points[pickedi].first, newq.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            E += V_ele;
        }

        E += 0.5*(qdot_1 - qdot).transpose()*M*(qdot_1 - qdot);

        return E;
    };

    auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) {

        assemble_forces(f, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C,D, energy_type);

        for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
            dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose()*q2+x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            f.segment<3>(3*constant_pickeds[pickedi]) -= dV_mouse.segment<3>(3);
        }

        f = P*f;
    };

    //assemble stiffness matrix,
    auto stiffness = [&](Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2) {
        assemble_stiffness(K, P.transpose()*q2+x0, P.transpose()*qdot2, V, T, v0, C, D, energy_type);
        K = P*K*P.transpose();
    };

    if(fully_implicit)
        implicit_euler(q, qdot, dt, M, energy, force, stiffness, tmp_qdot, tmp_force, tmp_stiffness);
    else
        linearly_implicit_euler(q, qdot, dt, M, force, stiffness, tmp_force, tmp_stiffness);



    KE = 0;
    PE = 0;
    VG = 0;

    for(unsigned int ei=0; ei<T.rows(); ++ei) {
        T_linear_tetrahedron(T_ele, P.transpose()*qdot, T.row(ei), density, v0(ei));
        KE += T_ele;

        V_linear_tetrahedron(V_ele, P.transpose()*q+x0, V, T.row(ei), v0(ei), C, D, energy_type);
        PE += V_ele;

        calc_volume(VG_ele,P.transpose()*q+x0, V, T.row(ei));
        VG += VG_ele;
    }
    if(initial < 10 or pow(10,2*oom_ke)< KE or pow(10,2*oom_pe)< PE){
        oom_ke = floor(log10(abs(KE)));
        oom_pe = floor(log10(abs(PE)));
        initial +=1;
    }
    string bun = "";
    if(bunny){
        bun = "bunny";
    }else{
        bun="arma";
    }
    KE = KE / ((range_ke_map[energy_type +"_"+bun])*fm);
    PE = PE / (range_pe_map[energy_type +"_"+bun]*fm);
    VG = VG * 1e2/ v0.sum();
    Visualize::add_energy_2(t, KE, PE, VG);
     std::cout << t << "\t"<< VG<<std::endl;
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    Visualize::update_vertex_positions(0, P.transpose()*q + x0);

}
std::string next_energy(){
    std::string now = energies[cur_energy];
    cur_energy += 1;
    return now;
}
bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    if(key == 'S') {

        skinning_on = !skinning_on;
        Visualize::toggle_skinning(skinning_on);
    }else if(key == 'E'){
        energy_type = next_energy();
    }

    return false;
}
void parse_options(int argc, char** argv){
    for(int i=1; i<argc-1 ;i+=2){
        if(strcmp(argv[i], "--mesh") == 0 || strcmp(argv[i], "-m") == 0){

            if(strcmp(argv[i+1], "bunny") == 0) {
                igl::readMESH("../data/coarser_bunny.mesh",V,T, F);
                igl::readOBJ("../data/bunny_skin.obj", V_skin, F_skin);
                bunny = true;
                fully_implicit = true;
            }else{
                read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");
                igl::readOBJ("../data/armadillo.obj", V_skin, F_skin);

                bunny = false;
                fully_implicit = true;

            }
        }else if(strcmp(argv[i], "--poisson_rate") == 0|| strcmp(argv[i], "-mu")== 0){
            if(atof(argv[i+1]) > 0.3 && atof(argv[i+1]) < 0.49) {
                mu = atof(argv[i+1]); //poissons ratio
                std::cout<<"Using poissons ratio: " << mu << "\n";
            }
            else{
                mu = 0.4; //poissons ratio
                std::cout<<"Using default poissons ratio: " << mu << "\n";
            }
        }else if(strcmp(argv[i], "--young_modulus")== 0 || strcmp(argv[i], "-ym")== 0){
            if(atof(argv[i+1]) > 10000.0 && atof(argv[i+1]) < 10000000.0) {
                YM = atof(argv[i+1]); //young's modulus
                std::cout<<"Using Young Modulus: " << YM << "\n";
            }
            else{
                YM = 6e5; //young's modulus
                std::cout<<"Using default Young Modulus: " << YM << "\n";
            }
        }else if(strcmp(argv[i], "--force_multiplier") == 0|| strcmp(argv[i], "-fm")== 0){
            if(atof(argv[i+1]) > 1.0 && atof(argv[i+1]) < 50.0) {
                fm = atof(argv[i+1]);
                std::cout<<"Using Force Multiplier: " << fm << "\n";
            }
            else{
                fm = 1.0;
                std::cout<<"Using default Force Multiplier: " << fm << "\n";
            }
        }else if(strcmp(argv[i], "--energy_model") == 0|| strcmp(argv[i], "-em")== 0){
            if(strcmp(argv[i+1],energies[0]) ==0 || strcmp(argv[i+1],energies[1])==0 || strcmp(argv[i+1],energies[2])==0 || strcmp(argv[i+1],energies[3])==0 || strcmp(argv[i+1],energies[4])==0){
                std::cout<<"Using Energy Model: " << argv[i+1] << "\n";
                energy_type = argv[i+1];
            }else{
                std::cout<<"Using default Energy Model: " << energies[0] << "\n";
                energy_type = energies[0];
            }
        }else if(strcmp(argv[i], "--sim_volume") == 0|| strcmp(argv[i], "-sv")== 0){
            if(strcmp(argv[i+1],"true") ==0 ){
                std::cout<<"Simulating Volume Gain " << argv[i+1] << "\n";
                simulate_vol = true;
            }



    }
}
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 
    read_tetgen(V,T, "../data/arma_6.node", "../data/arma_6.ele");
    igl::readOBJ("../data/armadillo.obj", V_skin, F_skin);

    bunny = false;
    fully_implicit = true;


    if(argc > 1) {
         parse_options(argc, argv);
    }
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();
    
    build_skinning_matrix(N, V, T, V_skin);

    //setup simulation 
    init_state(q,qdot,V);

    //add geometry to scene
    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::toggle_skinning(true);
    string bun = "";
    //bunny
    if(bunny){
        Visualize::set_picking_tolerance(1.);
        bun = "bunny";
    }else{
        Visualize::set_picking_tolerance(0.01);
        bun = "arma";
    }

    //volumes of all elements
    igl::volume(V,T, v0);

    //Mass Matrix
    mass_matrix_mesh(M, qdot, T, density, v0);
    
    if(M.rows() == 0) {
        std::cout<<"Mass Matrix not implemented, quitting \n";
        std::exit(0);
    }
    
    //setup constraint matrix
    if(bunny)
        find_min_vertices(fixed_point_indices, V, 3);
    else
        find_min_vertices(fixed_point_indices, V, 0.1);

    //material properties
    //bunny
    if(bunny) {
        YM = 6e6; //young's modulus
        mu = 0.4; //poissons ratio
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e8;
    } else {
        //arma
        D = 0.5*(YM*mu)/((1.0+mu)*(1.0-2.0*mu));
        C = 0.5*YM/(2.0*(1.0+mu));
        k_selected = 1e5;

    }

    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixed_point_constraints(P, q.rows(), fixed_point_indices);
    
    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else    
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();

    //igl additional menu setup
    // Add content to the default menu window

    Visualize::viewer_menu().callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * Visualize::viewer_menu().menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
            "Energy Plot", nullptr,
            ImGuiWindowFlags_NoSavedSettings

        );

        ImVec2 min = ImGui::GetWindowContentRegionMin();
        ImVec2 max = ImGui::GetWindowContentRegionMax();

        max.x = ( max.x - min.x ) / 2;
        max.y -= min.y + ImGui::GetItemsLineHeightWithSpacing() * 3;

        Visualize::plot_energy("T", 1, ImVec2(-15,10), ImVec2(0,10), ImGui::GetColorU32(ImGuiCol_PlotLines));
        Visualize::plot_energy("V", 2, ImVec2(-15,10), ImVec2(-30, 30), ImGui::GetColorU32(ImGuiCol_HeaderActive));

        if(simulate_vol){
            Visualize::plot_energy("Volume", 3, ImVec2(-15,10), ImVec2(80,120), ImGui::GetColorU32(ImGuiCol_ColumnActive));
        }else{
            Visualize::plot_energy("T+V", 3, ImVec2(-15,10), ImVec2(-30,30), ImGui::GetColorU32(ImGuiCol_ColumnActive));
        }

        ImGui::End();
    };

    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

