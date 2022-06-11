#include "Scancontext.h"

// namespace SC2
// {

void coreImportTest (void)
{
    cout << "scancontext lib is successfully imported." << endl;
} // coreImportTest


float rad2deg(float radians)
{
    return radians * 180.0 / M_PI;
}

float deg2rad(float degrees)
{
    return degrees * M_PI / 180.0;
}


float xy2theta( const float & _x, const float & _y )
{
    if ( (_x >= 0) & (_y >= 0)) 
        return (180/M_PI) * atan(_y / _x);

    if ( (_x < 0) & (_y >= 0)) 
        return 180 - ( (180/M_PI) * atan(_y / (-_x)) );

    if ( (_x < 0) & (_y < 0)) 
        return 180 + ( (180/M_PI) * atan(_y / _x) );

    if ( (_x >= 0) & (_y < 0))
        return 360 - ( (180/M_PI) * atan((-_y) / _x) );
} // xy2theta


/**
 * @brief 对SC进行旋转
 * 
 * @param _mat 待旋转的SC
 * @param _num_shift 旋转的sector数量
 * @return MatrixXd 旋转后的SC
 */
MatrixXd circshift( MatrixXd &_mat, int _num_shift )
{
    // shift columns to right direction 
    assert(_num_shift >= 0);

    if( _num_shift == 0 )
    {
        MatrixXd shifted_mat( _mat );
        return shifted_mat; // Early return 
    }

    MatrixXd shifted_mat = MatrixXd::Zero( _mat.rows(), _mat.cols() );
    for ( int col_idx = 0; col_idx < _mat.cols(); col_idx++ )
    {
        int new_location = (col_idx + _num_shift) % _mat.cols();
        shifted_mat.col(new_location) = _mat.col(col_idx);
    }

    return shifted_mat;

} // circshift


/**
 * @brief 把矩阵转成向量，eig2stdvec就是Eigen to std::vector的意思
 * 
 * @param _eigmat[in] 矩阵 
 * @return std::vector<float> 浮点向量
 */
std::vector<float> eig2stdvec( MatrixXd _eigmat )
{
    std::vector<float> vec( _eigmat.data(), _eigmat.data() + _eigmat.size() );
    return vec;
} // eig2stdvec


/**
 * @brief 计算两个SC的距离
 * 
 * 基本方法是sector的向量夹角的余弦值，向量夹角越小，相似度越高
 * 
 * 1. 同步遍历两个SC的每一个sector
 * 2. 忽略全零的sector
 * 3. 使用内积法计算两个sector向量夹角的余弦值
 * 4. 累加所有sector向量夹角的余弦值
 * 5. 求向量夹角余弦的均值，然后求相对于1的补，获得两个SC的距离
 * 
 * @param _sc1 
 * @param _sc2 
 * @return double 
 */
double SCManager::distDirectSC ( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    /* 同步遍历两个SC的每一个sector */
    int num_eff_cols = 0; // i.e., to exclude all-nonzero sector
    double sum_sector_similarity = 0;
    for ( int col_idx = 0; col_idx < _sc1.cols(); col_idx++ )
    {
        VectorXd col_sc1 = _sc1.col(col_idx);
        VectorXd col_sc2 = _sc2.col(col_idx);
        
        /* 忽略全零的sector */
        if( (col_sc1.norm() == 0) | (col_sc2.norm() == 0) )
            continue; // don't count this sector pair. 

        /* 使用内积法求两个sector向量夹角的cos值，相似度越高，越接近于1 */
        double sector_similarity = col_sc1.dot(col_sc2) / (col_sc1.norm() * col_sc2.norm());

        /* 累加所有sector向量夹角的余弦值 */
        sum_sector_similarity = sum_sector_similarity + sector_similarity;
        num_eff_cols = num_eff_cols + 1;
    }
    /* 计算sector向量夹角余弦的均值，取均值相对于1的补，即相似度越高，补越接近于0，越小 */
    double sc_sim = sum_sector_similarity / num_eff_cols;
    return 1.0 - sc_sim;

} // distDirectSC


/**
 * @brief 使用sector key快速对齐两个SC
 * 
 *  SC的旋转不变性就在这里体现了出来
 * 
 * 1. 对vkey2进行旋转，每旋转一格，求vkey1与旋转后vkey2的差值，旋转一圈
 * 2. 记录并返回最小的差值以及对应的vkey2的旋转index
 * 
 * NOTE: 对于环形雷达，可以使用旋转并求差值的方法对齐两个SC；
 *       对于面阵雷达，由于视角有限，旋转后两个SC的sector key将无法对齐。
 * 
 * @param _vkey1 SC1的sector key
 * @param _vkey2 SC2的sector key
 * @return int _vkey1的旋转index
 */
int SCManager::fastAlignUsingVkey( MatrixXd & _vkey1, MatrixXd & _vkey2)
{
    int argmin_vkey_shift = 0;
    double min_veky_diff_norm = 10000000;
    for ( int shift_idx = 0; shift_idx < _vkey1.cols(); shift_idx++ )
    {
        /* 旋转vkey2 */
        MatrixXd vkey2_shifted = circshift(_vkey2, shift_idx);

        /* vkey1减旋转后的vkey2, 求差*/
        MatrixXd vkey_diff = _vkey1 - vkey2_shifted;

        /* 记录最小的差值，以及旋转的index */
        double cur_diff_norm = vkey_diff.norm();
        if( cur_diff_norm < min_veky_diff_norm )
        {
            argmin_vkey_shift = shift_idx;
            min_veky_diff_norm = cur_diff_norm;
        }
    }

    /* 返回旋转后最小的差值，以及旋转的index */
    return argmin_vkey_shift;

} // fastAlignUsingVkey


/**
 * @brief 计算两个SC的距离
 * 
 * 1. 分别取得两个SC的sector key
 * 2. 通过旋转SC2的sector key并求差值的方法对齐两个SC，获得SC2的旋转index
 * 3. 将旋转index及其两侧半径为3的index，总共7个index添加到搜索空间
 * 4. 遍历搜索空间中的7个index
 *   a. 使用旋转index对SC2进行旋转
 *   b. 旋转后计算两个SC的距离，使用的是sector向量夹角的余弦值
 *   c. 找出SC距离的最小值
 * 5. 返回SC距离的最小值以及对应的旋转index
 * 
 * @param _sc1 
 * @param _sc2 
 * @return std::pair<double, int> 两个SC的距离
 */
std::pair<double, int> SCManager::distanceBtnScanContext( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    // 1. fast align using variant key (not in original IROS18)

    /* 分别取得两个SC的sector key */
    MatrixXd vkey_sc1 = makeSectorkeyFromScancontext( _sc1 );
    MatrixXd vkey_sc2 = makeSectorkeyFromScancontext( _sc2 );
    /* 利用SC的旋转不变性，对齐两个SC，获得sc2的旋转index */
    int argmin_vkey_shift = fastAlignUsingVkey( vkey_sc1, vkey_sc2 );

    /* 搜索半径 = 0.5 * 0.1 * 60 = 3 */
    const int SEARCH_RADIUS = round( 0.5 * SEARCH_RATIO * _sc1.cols() ); // a half of search range 
    /* 将SC1和SC2对齐的初始旋转index添加到搜索空间 */
    std::vector<int> shift_idx_search_space { argmin_vkey_shift };
    /* 将初始index及其两侧以3为半径的6个index添加到搜索空间 */
    for ( int ii = 1; ii < SEARCH_RADIUS + 1; ii++ )
    {
        shift_idx_search_space.push_back( (argmin_vkey_shift + ii + _sc1.cols()) % _sc1.cols() );
        shift_idx_search_space.push_back( (argmin_vkey_shift - ii + _sc1.cols()) % _sc1.cols() );
    }
    /* 对搜索空间进行排序，此时总共有7个旋转index */
    std::sort(shift_idx_search_space.begin(), shift_idx_search_space.end());

    // 2. fast columnwise diff 

    /* 遍历搜索空间中的7个旋转index */
    int argmin_shift = 0;
    double min_sc_dist = 10000000;
    for ( int num_shift: shift_idx_search_space )
    {
        /* 使用旋转index对SC2进行旋转 */
        MatrixXd sc2_shifted = circshift(_sc2, num_shift);

        /* 旋转后计算两个SC的距离，使用的是sector向量夹角的余弦值 */
        double cur_sc_dist = distDirectSC( _sc1, sc2_shifted );

        /* 找出SC距离的最小值 */
        if( cur_sc_dist < min_sc_dist )
        {
            argmin_shift = num_shift;
            min_sc_dist = cur_sc_dist;
        }
    }

    /* 返回最小距离值 */
    return make_pair(min_sc_dist, argmin_shift);

} // distanceBtnScanContext


/**
 * @brief 读取单帧点云数据，建立Scan-Context
 * 
 * 1. 建立由20个ring（环）和60个sector（扇区）组成的环形网格，即Scan-Context，每个网格记录最大高度值
 * 2. 遍历点云中所有的点
 *   a. 将点的坐标转成极坐标，计算出每个点对应的网格，点的半径超过80米的点丢弃
 *   b. 记录每个网格中z轴的最大值，即最高点
 * 3. 没有点的网格设置为0
 * 
 * @param _scan_down 关键帧点云
 * @return MatrixXd 返回建立好的Scan-Context
 */
MatrixXd SCManager::makeScancontext( pcl::PointCloud<SCPointType> & _scan_down )
{
    TicToc t_making_desc;

    int num_pts_scan_down = _scan_down.points.size();

    /* 建立由Ring和Sector组成的SC容器——矩阵 */
    // main
    const int NO_POINT = -1000;
    MatrixXd desc = NO_POINT * MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);

    /* 遍历点云中所有的点 */
    SCPointType pt;
    float azim_angle, azim_range; // wihtin 2d plane
    int ring_idx, sctor_idx;
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z = _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // naive adding is ok (all points should be > 0).

        /* 笛卡尔坐标转极坐标 */
        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);

        /* 点的半径超过80米，丢弃 */
        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )
            continue;

        /* 根据轴长计算出ring index */
        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
        /* 根据角度计算出sector index */
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );

        /* 记录每个网格的最大z值 */
        // taking maximum z 
        if ( desc(ring_idx-1, sctor_idx-1) < pt.z ) // -1 means cpp starts from 0
            desc(ring_idx-1, sctor_idx-1) = pt.z; // update for taking maximum value at that bin
    }

    /* 没有点的网格设置为0 */
    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )
            if( desc(row_idx, col_idx) == NO_POINT )
                desc(row_idx, col_idx) = 0;

    t_making_desc.toc("PolarContext making");

    return desc;
} // SCManager::makeScancontext


/**
 * @brief 计算出当前关键帧点云的ring key
 * 
 * 1. 以Scan-Context为输入，计算出每一行（环）的高度均值
 * 2. 将20个环的高度均值组成向量，作为当前帧点云的ring key
 * 
 * @param _desc 当前帧点云对应的Scan-Context
 * @return MatrixXd 当前帧点云对应的ring key
 */
MatrixXd SCManager::makeRingkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: rowwise mean vector
    */
    Eigen::MatrixXd invariant_key(_desc.rows(), 1);
    for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ )
    {
        Eigen::MatrixXd curr_row = _desc.row(row_idx);
        invariant_key(row_idx, 0) = curr_row.mean();
    }

    return invariant_key;
} // SCManager::makeRingkeyFromScancontext


/**
 * @brief 计算出当前关键帧点云的sector key
 * 
 * 1. 以Scan-Context为输入，计算出每一列（扇区）的高度均值
 * 2. 将60个扇区的高度均值组成向量，作为当前帧点云的sector key
 * 
 * @param _desc 当前帧点云对应的Scan-Context
 * @return MatrixXd 当前帧点云对应的sector key
 */
MatrixXd SCManager::makeSectorkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: columnwise mean vector
    */
    Eigen::MatrixXd variant_key(1, _desc.cols());
    for ( int col_idx = 0; col_idx < _desc.cols(); col_idx++ )
    {
        Eigen::MatrixXd curr_col = _desc.col(col_idx);
        variant_key(0, col_idx) = curr_col.mean();
    }

    return variant_key;
} // SCManager::makeSectorkeyFromScancontext


const Eigen::MatrixXd& SCManager::getConstRefRecentSCD(void)
{
    return polarcontexts_.back();
}

/**
 * @brief 将关键帧点云添加到SC探测器
 * 
 * 1. 计算出当前关键帧点云的Scan-Context
 * 2. 计算出当前关键帧点云的ring key
 * 3. 计算出当前关键帧点云的sector key
 * 4. 将ring key转成浮点向量格式ring key vec
 * 
 * @param _scan_down 关键帧点云
 */
void SCManager::makeAndSaveScancontextAndKeys( pcl::PointCloud<SCPointType> & _scan_down )
{
    /* 分别计算当前关键帧的Scan-Context、ringkey、sectorkey等，sectorkey没有用 */
    Eigen::MatrixXd sc = makeScancontext(_scan_down); // v1 
    Eigen::MatrixXd ringkey = makeRingkeyFromScancontext( sc );
    Eigen::MatrixXd sectorkey = makeSectorkeyFromScancontext( sc );
    /* 将ringkey转成浮点向量的形式 */
    std::vector<float> polarcontext_invkey_vec = eig2stdvec( ringkey );

    /* 保存sc、ringkey、sectorkey和ringkey vec，ringkey和sectorkey都没有用，只用ringkey vec */
    polarcontexts_.push_back( sc ); 
    polarcontext_invkeys_.push_back( ringkey );
    polarcontext_vkeys_.push_back( sectorkey );
    polarcontext_invkeys_mat_.push_back( polarcontext_invkey_vec );

    // cout <<polarcontext_vkeys_.size() << endl;

} // SCManager::makeAndSaveScancontextAndKeys


/**
 * @brief 检测可能存在的回环
 * 
 * 总的思路是用ringkey的旋转不变性找出候选关键帧，然后用sectorkey确定回环帧，以及当前帧与回环帧的航向差值（yaw）。
 * 
 * 1. 取得当前关键帧的ringkey向量和SC
 * 2. 为所有历史关键帧的ringkey向量建立高维空间kdtree，每个关键帧的ringkey向量都是kdtree中的一个点
 * 3. 在kdtree中搜索与当前关键帧的ringkey向量相似的候选关键帧
 * 4. 遍历候选关键帧，计算候选关键帧与当前关键帧SC的距离，找出距离最小的候选关键帧，返回候选关键帧的id
 * 
 * @return std::pair<int, float> 
 */
std::pair<int, float> SCManager::detectLoopClosureID ( void )
{
    int loop_id { -1 }; // init with -1, -1 means no loop (== LeGO-LOAM's variable "closestHistoryFrameID")

    /* 取得当前关键帧对应的SC和ringkey_vec */
    auto curr_key = polarcontext_invkeys_mat_.back(); // current observation (query)
    auto curr_desc = polarcontexts_.back(); // current observation (query)

    /* 
     * step 1: candidates from ringkey tree_
     */

    /* 总的关键帧数量不足30个，直接返回 */
    if( (int)polarcontext_invkeys_mat_.size() < NUM_EXCLUDE_RECENT + 1)
    {
        std::pair<int, float> result {loop_id, 0.0};
        return result; // Early return 
    }

    /* 用既往所有关键帧的ringkey_vec构造kdtree，不需要周期都进行，每隔10帧进行一次 */
    // tree_ reconstruction (not mandatory to make everytime)
    if( tree_making_period_conter % TREE_MAKING_PERIOD_ == 0) // to save computation cost
    {
        TicToc t_tree_construction;

        /* 把历史上的ringkey_vec都取出来，放到polarcontext_invkeys_to_search_中 */
        polarcontext_invkeys_to_search_.clear();
        polarcontext_invkeys_to_search_.assign( polarcontext_invkeys_mat_.begin(), polarcontext_invkeys_mat_.end() - NUM_EXCLUDE_RECENT ) ;

        /**
         * 建立kdtree
         * 
         * 这里的kdtree不是通常的三维空间点，而是高维空间点，每一个关键帧对应的ringkey就是高维空间中的一个点，
         * 维度是ringkey向量的长度20。
         * 
         * 由于ringkey向量由关键帧的多个ring的均值组成，因此具有旋转不变性，即只要两个关键帧在同一个位置，
         * 不论方向如何，ringkey向量的值应该都是一样的。
         */
        polarcontext_tree_.reset(); 
        polarcontext_tree_ = std::make_unique<InvKeyTree>(PC_NUM_RING /* dim */, polarcontext_invkeys_to_search_, 10 /* max leaf */ );
        // tree_ptr_->index->buildIndex(); // inernally called in the constructor of InvKeyTree (for detail, refer the nanoflann and KDtreeVectorOfVectorsAdaptor)
        t_tree_construction.toc("Tree construction");
    }
    tree_making_period_conter = tree_making_period_conter + 1;
        
    double min_dist = 10000000; // init with somthing large
    int nn_align = 0;
    int nn_idx = 0;

    // knn search
    std::vector<size_t> candidate_indexes( NUM_CANDIDATES_FROM_TREE ); 
    std::vector<float> out_dists_sqr( NUM_CANDIDATES_FROM_TREE );

    /* 在kdtree中搜索最接近的（即最相似）的ringkey，获得候选帧的index */
    TicToc t_tree_search;
    nanoflann::KNNResultSet<float> knnsearch_result( NUM_CANDIDATES_FROM_TREE );
    knnsearch_result.init( &candidate_indexes[0], &out_dists_sqr[0] );
    polarcontext_tree_->index->findNeighbors( knnsearch_result, &curr_key[0] /* query */, nanoflann::SearchParams(10) ); 
    t_tree_search.toc("Tree search");

    /* 
     *  step 2: pairwise distance (find optimal columnwise best-fit using cosine distance)
     */

    /* 遍历候选帧，计算候选帧的SC与当前帧SC的空间距离，找出距离最小的候选帧 */    
    TicToc t_calc_dist;   
    for ( int candidate_iter_idx = 0; candidate_iter_idx < NUM_CANDIDATES_FROM_TREE; candidate_iter_idx++ )
    {
        /* 通过候选者index取得对应的SC */
        MatrixXd polarcontext_candidate = polarcontexts_[ candidate_indexes[candidate_iter_idx] ];
        /* 计算当前帧SC和候选帧SC的距离 */
        std::pair<double, int> sc_dist_result = distanceBtnScanContext( curr_desc, polarcontext_candidate ); 
        
        double candidate_dist = sc_dist_result.first;
        int candidate_align = sc_dist_result.second;

        /* 记录距离最小的候选帧 */
        if( candidate_dist < min_dist )
        {
            min_dist = candidate_dist;
            nn_align = candidate_align;

            nn_idx = candidate_indexes[candidate_iter_idx];
        }
    }
    t_calc_dist.toc("Distance calc");

    /* 如果SC的最小距离低于阈值，则返回结果，否则返回-1 */
    /* 
     * loop threshold check
     */
    if( min_dist < SC_DIST_THRES )
    {
        loop_id = nn_idx; 
    
        // std::cout.precision(3); 
        cout << "[Loop found] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        cout << "[Loop found] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }
    else
    {
        std::cout.precision(3); 
        cout << "[Not loop] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        cout << "[Not loop] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }

    // To do: return also nn_align (i.e., yaw diff)
    float yaw_diff_rad = deg2rad(nn_align * PC_UNIT_SECTORANGLE);
    std::pair<int, float> result {loop_id, yaw_diff_rad};

    return result;

} // SCManager::detectLoopClosureID

// } // namespace SC2