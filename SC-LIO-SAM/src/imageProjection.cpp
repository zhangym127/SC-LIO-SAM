#include "utility.h"
#include "lio_sam/cloud_info.h"

struct VelodynePointXYZIRT
{
    PCL_ADD_POINT4D
    PCL_ADD_INTENSITY;
    uint16_t ring;
    float time;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;
POINT_CLOUD_REGISTER_POINT_STRUCT (VelodynePointXYZIRT,
    (float, x, x) (float, y, y) (float, z, z) (float, intensity, intensity)
    (uint16_t, ring, ring) (float, time, time)
)

struct OusterPointXYZIRT {
    PCL_ADD_POINT4D;
    float intensity;
    uint32_t t;
    uint16_t reflectivity;
    uint8_t ring;
    uint16_t noise;
    uint32_t range;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;
POINT_CLOUD_REGISTER_POINT_STRUCT(OusterPointXYZIRT,
    (float, x, x) (float, y, y) (float, z, z) (float, intensity, intensity)
    (uint32_t, t, t) (uint16_t, reflectivity, reflectivity)
    (uint8_t, ring, ring) (uint16_t, noise, noise) (uint32_t, range, range)
)

struct MulranPointXYZIRT { // from the file player's topic https://github.com/irapkaist/file_player_mulran, see https://github.com/irapkaist/file_player_mulran/blob/17da0cb6ef66b4971ec943ab8d234aa25da33e7e/src/ROSThread.cpp#L7
     PCL_ADD_POINT4D;
     float intensity;
     uint32_t t;
     int ring;

     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 }EIGEN_ALIGN16;
 POINT_CLOUD_REGISTER_POINT_STRUCT (MulranPointXYZIRT,
     (float, x, x) (float, y, y) (float, z, z) (float, intensity, intensity)
     (uint32_t, t, t) (int, ring, ring)
 )

// Use the Velodyne point format as a common representation
using PointXYZIRT = VelodynePointXYZIRT;

const int queueLength = 2000;

class ImageProjection : public ParamServer
{
private:

    std::mutex imuLock;
    std::mutex odoLock;

    ros::Subscriber subLaserCloud;
    ros::Publisher  pubLaserCloud;
    
    ros::Publisher pubExtractedCloud;
    ros::Publisher pubLaserCloudInfo;

    ros::Subscriber subImu;
    std::deque<sensor_msgs::Imu> imuQueue;

    ros::Subscriber subOdom;
    std::deque<nav_msgs::Odometry> odomQueue;

    std::deque<sensor_msgs::PointCloud2> cloudQueue;
    sensor_msgs::PointCloud2 currentCloudMsg;

    double *imuTime = new double[queueLength];
    double *imuRotX = new double[queueLength];
    double *imuRotY = new double[queueLength];
    double *imuRotZ = new double[queueLength];

    int imuPointerCur;
    bool firstPointFlag;
    Eigen::Affine3f transStartInverse;

    pcl::PointCloud<PointXYZIRT>::Ptr laserCloudIn;
    pcl::PointCloud<OusterPointXYZIRT>::Ptr tmpOusterCloudIn;
    pcl::PointCloud<MulranPointXYZIRT>::Ptr tmpMulranCloudIn;
    pcl::PointCloud<PointType>::Ptr   fullCloud;
    pcl::PointCloud<PointType>::Ptr   extractedCloud;

    int deskewFlag;
    cv::Mat rangeMat;

    bool odomDeskewFlag;
    float odomIncreX;
    float odomIncreY;
    float odomIncreZ;

    lio_sam::cloud_info cloudInfo;
    double timeScanCur;
    double timeScanEnd;
    std_msgs::Header cloudHeader;


public:

	/**
	 * 构造函数 
	 */
    ImageProjection():
    deskewFlag(0)
    {
		/* 订阅IMU原始数据，交给imuHandler处理 */
        subImu        = nh.subscribe<sensor_msgs::Imu>(imuTopic, 2000, &ImageProjection::imuHandler, this, ros::TransportHints().tcpNoDelay());
		/* 订阅IMUPreintegration发布的里程计信息，交给odometryHandler处理 */
        subOdom       = nh.subscribe<nav_msgs::Odometry>(odomTopic+"_incremental", 2000, &ImageProjection::odometryHandler, this, ros::TransportHints().tcpNoDelay());
        /* 订阅点云的原始数据，交给cloudHandler处理，cloudHandler是本类的核心函数 */
        subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2>(pointCloudTopic, 5, &ImageProjection::cloudHandler, this, ros::TransportHints().tcpNoDelay());

		/* 创建发布，准备发布点云数据 */
        pubExtractedCloud = nh.advertise<sensor_msgs::PointCloud2> ("lio_sam/deskew/cloud_deskewed", 1);
        pubLaserCloudInfo = nh.advertise<lio_sam::cloud_info> ("lio_sam/deskew/cloud_info", 1);

        allocateMemory();
        resetParameters();

        pcl::console::setVerbosityLevel(pcl::console::L_ERROR);
    }

    void allocateMemory()
    {
        laserCloudIn.reset(new pcl::PointCloud<PointXYZIRT>());
        tmpOusterCloudIn.reset(new pcl::PointCloud<OusterPointXYZIRT>());
        tmpMulranCloudIn.reset(new pcl::PointCloud<MulranPointXYZIRT>());
        fullCloud.reset(new pcl::PointCloud<PointType>());
        extractedCloud.reset(new pcl::PointCloud<PointType>());

        fullCloud->points.resize(N_SCAN*Horizon_SCAN);

        cloudInfo.startRingIndex.assign(N_SCAN, 0);
        cloudInfo.endRingIndex.assign(N_SCAN, 0);

        cloudInfo.pointColInd.assign(N_SCAN*Horizon_SCAN, 0);
        cloudInfo.pointRange.assign(N_SCAN*Horizon_SCAN, 0);

        resetParameters();
    }

    void resetParameters()
    {
        laserCloudIn->clear();
        extractedCloud->clear();
        // reset range matrix for range image projection
        rangeMat = cv::Mat(N_SCAN, Horizon_SCAN, CV_32F, cv::Scalar::all(FLT_MAX));

        imuPointerCur = 0;
        firstPointFlag = true;
        odomDeskewFlag = false;

        for (int i = 0; i < queueLength; ++i)
        {
            imuTime[i] = 0;
            imuRotX[i] = 0;
            imuRotY[i] = 0;
            imuRotZ[i] = 0;
        }
    }

    ~ImageProjection(){}

    /**
     * @brief 获得IMU的原始数据，压入imuQueue队列
     * 
     * 注意：这里取得的IMU原始数据与IMU预积分取得的完全相同
     * 
     * @param imuMsg 一帧IMU原始数据
     */
    void imuHandler(const sensor_msgs::Imu::ConstPtr& imuMsg)
    {
        sensor_msgs::Imu thisImu = imuConverter(*imuMsg);

        std::lock_guard<std::mutex> lock1(imuLock);
        imuQueue.push_back(thisImu);

        // debug IMU data
        // cout << std::setprecision(6);
        // cout << "IMU acc: " << endl;
        // cout << "x: " << thisImu.linear_acceleration.x << 
        //       ", y: " << thisImu.linear_acceleration.y << 
        //       ", z: " << thisImu.linear_acceleration.z << endl;
        // cout << "IMU gyro: " << endl;
        // cout << "x: " << thisImu.angular_velocity.x << 
        //       ", y: " << thisImu.angular_velocity.y << 
        //       ", z: " << thisImu.angular_velocity.z << endl;
        // double imuRoll, imuPitch, imuYaw;
        // tf::Quaternion orientation;
        // tf::quaternionMsgToTF(thisImu.orientation, orientation);
        // tf::Matrix3x3(orientation).getRPY(imuRoll, imuPitch, imuYaw);
        // cout << "IMU roll pitch yaw: " << endl;
        // cout << "roll: " << imuRoll << ", pitch: " << imuPitch << ", yaw: " << imuYaw << endl << endl;
    }

	/**
	 * 处理来自IMUPreintegration的里程计数据，仅仅是压入odomQueue
	 */
    void odometryHandler(const nav_msgs::Odometry::ConstPtr& odometryMsg)
    {
        std::lock_guard<std::mutex> lock2(odoLock);
        odomQueue.push_back(*odometryMsg);
    }

    /**
     * @brief 处理点云的原始数据，对点云进行帧内校正，并按照扫描线号进行重排列
     * 
     * 本函数是imageProjection类的核心函数，也是原始点云订阅的处理函数
     * 
     * @param laserCloudMsg 
     */
    void cloudHandler(const sensor_msgs::PointCloud2ConstPtr& laserCloudMsg)
    {
		/* 缓存点云原始数据，并将点云从ROS转成PCL格式 */
        if (!cachePointCloud(laserCloudMsg))
            return;

		/* 获得当前点云帧对应的旋转和平移量，准备进行点云的帧内校正 */ 
        if (!deskewInfo())
            return;

		/* 按照行列对齐的方式对点云进行重排列，并进行帧内校正 */
        projectPointCloud();

		/* 从整理好的点云中提取有效点，一维化，准备发布 */
        cloudExtraction();

		/* 发布点云 */
        publishClouds();

        resetParameters();
    }

	/**
 	 * 缓存点云原始数据，并将点云从ROS转成PCL格式。
	 */
    bool cachePointCloud(const sensor_msgs::PointCloud2ConstPtr& laserCloudMsg)
    {
        // cache point cloud
        cloudQueue.push_back(*laserCloudMsg);
        if (cloudQueue.size() <= 2)
            return false;

        // convert cloud
        currentCloudMsg = std::move(cloudQueue.front());

        cloudHeader = currentCloudMsg.header;
		/* 当前帧点云TS和下一帧点云TS */
        timeScanCur = cloudHeader.stamp.toSec();

        cloudQueue.pop_front();

        if (sensor == SensorType::VELODYNE)
        {
            pcl::moveFromROSMsg(currentCloudMsg, *laserCloudIn);
        }
        else if (sensor == SensorType::OUSTER)
        {
            // Convert to Velodyne format
            pcl::moveFromROSMsg(currentCloudMsg, *tmpOusterCloudIn);
            laserCloudIn->points.resize(tmpOusterCloudIn->size());
            laserCloudIn->is_dense = tmpOusterCloudIn->is_dense;
            for (size_t i = 0; i < tmpOusterCloudIn->size(); i++)
            {
                auto &src = tmpOusterCloudIn->points[i];
                auto &dst = laserCloudIn->points[i];
                dst.x = src.x;
                dst.y = src.y;
                dst.z = src.z;
                dst.intensity = src.intensity;
                dst.ring = src.ring;
                dst.time = src.t * 1e-9f;
            }
        }
        else if (sensor == SensorType::MULRAN)
        {
            // Convert to Velodyne format
            pcl::moveFromROSMsg(currentCloudMsg, *tmpMulranCloudIn);
            laserCloudIn->points.resize(tmpMulranCloudIn->size());
            laserCloudIn->is_dense = tmpMulranCloudIn->is_dense;
            for (size_t i = 0; i < tmpMulranCloudIn->size(); i++)
            {
                auto &src = tmpMulranCloudIn->points[i];
                auto &dst = laserCloudIn->points[i];
                dst.x = src.x;
                dst.y = src.y;
                dst.z = src.z;
                dst.intensity = src.intensity;
                dst.ring = src.ring;
                dst.time = float(src.t);
            }
        }
        else
        {
            ROS_ERROR_STREAM("Unknown sensor type: " << int(sensor));
            ros::shutdown();
        }

        // get timestamp
        timeScanEnd = timeScanCur + laserCloudIn->points.back().time;

        // check dense flag
        if (laserCloudIn->is_dense == false)
        {
            ROS_ERROR("Point cloud is not in dense format, please remove NaN points first!");
            ros::shutdown();
        }

        // check ring channel
        static int ringFlag = 0;
        if (ringFlag == 0)
        {
            ringFlag = -1;
            for (int i = 0; i < (int)currentCloudMsg.fields.size(); ++i)
            {
                if (currentCloudMsg.fields[i].name == "ring")
                {
                    ringFlag = 1;
                    break;
                }
            }
            if (ringFlag == -1)
            {
                ROS_ERROR("Point cloud ring channel not available, please configure your point cloud data!");
                ros::shutdown();
            }
        }

        // check point time
        if (deskewFlag == 0)
        {
            deskewFlag = -1;
            for (auto &field : currentCloudMsg.fields)
            {
                if (field.name == "time" || field.name == "t")
                {
                    deskewFlag = 1;
                    break;
                }
            }
            if (deskewFlag == -1)
                ROS_WARN("Point cloud timestamp not available, deskew function disabled, system will drift significantly!");
        }

        return true;
    }

    /**
     * @brief 取得用于当前帧点云帧内校正的旋转和平移增量
     * 
     * 获得当前点云帧对应的旋转和平移量，旋转来自对IMU数据的积分，平移来自
	 * IMUPreintegration的里程计信息，用于点云的帧内校正。
     * 
     * @return true 
     * @return false 
     */
    bool deskewInfo()
    {
        std::lock_guard<std::mutex> lock1(imuLock);
        std::lock_guard<std::mutex> lock2(odoLock);

        // make sure IMU data available for the scan
        if (imuQueue.empty() || imuQueue.front().header.stamp.toSec() > timeScanCur || imuQueue.back().header.stamp.toSec() < timeScanEnd)
        {
            ROS_DEBUG("Waiting for IMU data ...");
            return false;
        }

        /* 取得当前点云帧对应的旋转，通过对IMU的角速度积分获得 */
        imuDeskewInfo();

        /* 取得当前点云帧对应的位移，从来自IMUPreintegration的里程计信息获得当前点云帧期间的总平移
         * 在当前点云中填入位姿的估计值，用于mapOptimization的scan to map优化 */
        odomDeskewInfo();

        return true;
    }

    /**
     * @brief 对当前点云帧对应的IMU角速度进行积分
     * 
     * 对IMU角速度进行积分，获得每一帧IMU数据相对于当前点云帧起始时刻对应的旋转量，
	 * 保存在imuRotX，imuRotY，imuRotZ中，用于点云的帧内校正。这个做法比较简单，
     * 与livox_horizon_loam使用李代数对IMU积分结果转成李代数进行插值的做法相比精
     * 度肯定会差一些，但是计算更简单。
     */
    void imuDeskewInfo()
    {
        cloudInfo.imuAvailable = false;

        while (!imuQueue.empty())
        {
            if (imuQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                imuQueue.pop_front();
            else
                break;
        }

        if (imuQueue.empty())
            return;

        imuPointerCur = 0;

        /* 遍历imuQueue中的每一帧IMU数据 */
        for (int i = 0; i < (int)imuQueue.size(); ++i)
        {
            sensor_msgs::Imu thisImuMsg = imuQueue[i];
            double currentImuTime = thisImuMsg.header.stamp.toSec();

            /** 
             * 取得当前帧IMU数据对应的roll，pitch和yaw
             * FIXME: 这里每取到一帧IMU数据都会去更新一遍cloudInfo中的roll，pitch和yaw不知道是啥意思
             * 难道IMU中的orientation四元数是积分后的当前帧位姿不成？
             */
            // get roll, pitch, and yaw estimation for this scan
            if (currentImuTime <= timeScanCur)
                imuRPY2rosRPY(&thisImuMsg, &cloudInfo.imuRollInit, &cloudInfo.imuPitchInit, &cloudInfo.imuYawInit);

            if (currentImuTime > timeScanEnd + 0.01)
                break;

            if (imuPointerCur == 0){
                imuRotX[0] = 0;
                imuRotY[0] = 0;
                imuRotZ[0] = 0;
                imuTime[0] = currentImuTime;
                ++imuPointerCur;
                continue;
            }

            /* 取得当前帧IMU数据的角速度 */
            // get angular velocity
            double angular_x, angular_y, angular_z;
            imuAngular2rosAngular(&thisImuMsg, &angular_x, &angular_y, &angular_z);

            /* 对角速度进行积分，获得当前帧IMU数据对应的旋转 */
            // integrate rotation
            double timeDiff = currentImuTime - imuTime[imuPointerCur-1];
            imuRotX[imuPointerCur] = imuRotX[imuPointerCur-1] + angular_x * timeDiff;
            imuRotY[imuPointerCur] = imuRotY[imuPointerCur-1] + angular_y * timeDiff;
            imuRotZ[imuPointerCur] = imuRotZ[imuPointerCur-1] + angular_z * timeDiff;
            imuTime[imuPointerCur] = currentImuTime;
            ++imuPointerCur;
        }

        --imuPointerCur;

        if (imuPointerCur <= 0)
            return;

        cloudInfo.imuAvailable = true;
    }

    /**
     * @brief 从来自IMUPreintegration的里程计信息获得当前点云帧期间的位移增量，用
	 * 于点云的帧内校正，保存在odomIncreX，odomIncreY，odomIncreZ三个变量中。
     * 在当前点云中填入位姿的估计值，用于mapOptimization的scan to map优化。
     */
    void odomDeskewInfo()
    {
        cloudInfo.odomAvailable = false;

        while (!odomQueue.empty())
        {
            if (odomQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                odomQueue.pop_front();
            else
                break;
        }

        if (odomQueue.empty())
            return;

        if (odomQueue.front().header.stamp.toSec() > timeScanCur)
            return;

        // get start odometry at the beinning of the scan
        nav_msgs::Odometry startOdomMsg;

        for (int i = 0; i < (int)odomQueue.size(); ++i)
        {
            startOdomMsg = odomQueue[i];

            if (ROS_TIME(&startOdomMsg) < timeScanCur)
                continue;
            else
                break;
        }

        tf::Quaternion orientation;
        tf::quaternionMsgToTF(startOdomMsg.pose.pose.orientation, orientation);

        double roll, pitch, yaw;
        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

        /* 在当前点云中填入位姿的估计值，用于mapOptimization的scan to map优化 */
        // Initial guess used in mapOptimization
        cloudInfo.initialGuessX = startOdomMsg.pose.pose.position.x;
        cloudInfo.initialGuessY = startOdomMsg.pose.pose.position.y;
        cloudInfo.initialGuessZ = startOdomMsg.pose.pose.position.z;
        cloudInfo.initialGuessRoll  = roll;
        cloudInfo.initialGuessPitch = pitch;
        cloudInfo.initialGuessYaw   = yaw;

        cloudInfo.odomAvailable = true;

        // get end odometry at the end of the scan
        odomDeskewFlag = false;

        if (odomQueue.back().header.stamp.toSec() < timeScanEnd)
            return;

        nav_msgs::Odometry endOdomMsg;

        for (int i = 0; i < (int)odomQueue.size(); ++i)
        {
            endOdomMsg = odomQueue[i];

            if (ROS_TIME(&endOdomMsg) < timeScanEnd)
                continue;
            else
                break;
        }

        if (int(round(startOdomMsg.pose.covariance[0])) != int(round(endOdomMsg.pose.covariance[0])))
            return;

        Eigen::Affine3f transBegin = pcl::getTransformation(startOdomMsg.pose.pose.position.x, startOdomMsg.pose.pose.position.y, startOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        tf::quaternionMsgToTF(endOdomMsg.pose.pose.orientation, orientation);
        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
        Eigen::Affine3f transEnd = pcl::getTransformation(endOdomMsg.pose.pose.position.x, endOdomMsg.pose.pose.position.y, endOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        Eigen::Affine3f transBt = transBegin.inverse() * transEnd;

        float rollIncre, pitchIncre, yawIncre;
        pcl::getTranslationAndEulerAngles(transBt, odomIncreX, odomIncreY, odomIncreZ, rollIncre, pitchIncre, yawIncre);

        odomDeskewFlag = true;
    }

    /**
     * @brief 获得指定时刻对应的旋转，旋转量来自IMU角速度积分imuRotX，imuRotY，imuRotZ
     * 
     * @param pointTime 
     * @param rotXCur 
     * @param rotYCur 
     * @param rotZCur 
     */
    void findRotation(double pointTime, float *rotXCur, float *rotYCur, float *rotZCur)
    {
        *rotXCur = 0; *rotYCur = 0; *rotZCur = 0;

        int imuPointerFront = 0;
        while (imuPointerFront < imuPointerCur)
        {
            if (pointTime < imuTime[imuPointerFront])
                break;
            ++imuPointerFront;
        }

        if (pointTime > imuTime[imuPointerFront] || imuPointerFront == 0)
        {
            *rotXCur = imuRotX[imuPointerFront];
            *rotYCur = imuRotY[imuPointerFront];
            *rotZCur = imuRotZ[imuPointerFront];
        } else {
            int imuPointerBack = imuPointerFront - 1;
            double ratioFront = (pointTime - imuTime[imuPointerBack]) / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
            double ratioBack = (imuTime[imuPointerFront] - pointTime) / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
            *rotXCur = imuRotX[imuPointerFront] * ratioFront + imuRotX[imuPointerBack] * ratioBack;
            *rotYCur = imuRotY[imuPointerFront] * ratioFront + imuRotY[imuPointerBack] * ratioBack;
            *rotZCur = imuRotZ[imuPointerFront] * ratioFront + imuRotZ[imuPointerBack] * ratioBack;
        }
    }

    /**
     * @brief 获得指定时刻对应的平移，平移量原本来自里程计，但是并未实施
     * 
     * FIXME: 原本该函数应该用odomIncreX，odomIncreY，odomIncreZ三个变量计算并输出当前点云帧对应的平移增量，
     * 但是全部被注释掉了，也就是点云帧校正的时候平移全部为零。
     * 
     * @param relTime 
     * @param posXCur 
     * @param posYCur 
     * @param posZCur 
     */
    void findPosition(double relTime, float *posXCur, float *posYCur, float *posZCur)
    {
        *posXCur = 0; *posYCur = 0; *posZCur = 0;

        // If the sensor moves relatively slow, like walking speed, positional deskew seems to have little benefits. Thus code below is commented.

        // if (cloudInfo.odomAvailable == false || odomDeskewFlag == false)
        //     return;

        // float ratio = relTime / (timeScanEnd - timeScanCur);

        // *posXCur = ratio * odomIncreX;
        // *posYCur = ratio * odomIncreY;
        // *posZCur = ratio * odomIncreZ;
    }

    /**
     * @brief 对当前点进行帧内校正，只校正了旋转，未校正平移
     * 
     * @param point[in] 待校正的当前点
     * @param relTime[in] 当前点对应的时间戳
     * @return PointType 返回校正的后点 
     */
    PointType deskewPoint(PointType *point, double relTime)
    {
        if (deskewFlag == -1 || cloudInfo.imuAvailable == false)
            return *point;

		/* 获得当前点的时间戳 */
        double pointTime = timeScanCur + relTime;

		/* 找到当前点对应的旋转，从imuRotX、imuRotY、imuRotZ中获得 */
        float rotXCur, rotYCur, rotZCur;
        findRotation(pointTime, &rotXCur, &rotYCur, &rotZCur);

		/* 找到当前点对应的平移，实际获取的平移是0 */
        float posXCur, posYCur, posZCur;
        findPosition(relTime, &posXCur, &posYCur, &posZCur);

		/* 获得当前帧第一个点位姿的逆 */
        if (firstPointFlag == true)
        {
            transStartInverse = (pcl::getTransformation(posXCur, posYCur, posZCur, rotXCur, rotYCur, rotZCur)).inverse();
            firstPointFlag = false;
        }

		/* 获得相对于当前帧第一个点位姿的位姿变换 */
        // transform points to start
        Eigen::Affine3f transFinal = pcl::getTransformation(posXCur, posYCur, posZCur, rotXCur, rotYCur, rotZCur);
        Eigen::Affine3f transBt = transStartInverse * transFinal;

		/* 对当前点进行位姿变换 */
        PointType newPoint;
        newPoint.x = transBt(0,0) * point->x + transBt(0,1) * point->y + transBt(0,2) * point->z + transBt(0,3);
        newPoint.y = transBt(1,0) * point->x + transBt(1,1) * point->y + transBt(1,2) * point->z + transBt(1,3);
        newPoint.z = transBt(2,0) * point->x + transBt(2,1) * point->y + transBt(2,2) * point->z + transBt(2,3);
        newPoint.intensity = point->intensity;

        return newPoint;
    }


    /**
     * @brief 对点云进行帧内校正，并投影到深度图，按照图像的扫描顺序进行重排序
     * 
     * 0. 从laserCloudIn中取到当前帧点云并遍历每一个点
     * 1. 根据点的扫描线号确定深度图的纵坐标
     * 2. 根据点的水平视角确定深度图的横坐标
     * 3. 将当前点映射到深度图，如果有重复点则剔除
     * 4. 对当前点进行帧内校正
     * 5. 更新当前点在深度图内的深度
     * 6. 根据深度图的顺序，对点云进行重排序，存储到fullCloud中
     * 
     * 这个代码是为Velodyne、Ouster等旋转结构的雷达打造的，特点在下面的代码中能够反映出来
     */
    void projectPointCloud()
    {
        int cloudSize = laserCloudIn->points.size();
        // range image projection
        for (int i = 0; i < cloudSize; ++i)
        {
            PointType thisPoint;
            thisPoint.x = laserCloudIn->points[i].x;
            thisPoint.y = laserCloudIn->points[i].y;
            thisPoint.z = laserCloudIn->points[i].z;
            thisPoint.intensity = laserCloudIn->points[i].intensity;

            /* 取得点到雷达坐标原点的距离，即深度 */
            float range = pointDistance(thisPoint);
            if (range < lidarMinRange || range > lidarMaxRange)
                continue;

			/* 获得当前点的行号（即激光雷达的扫描线号）对应到深度图的纵坐标rowIdn */
            int rowIdn = laserCloudIn->points[i].ring;
            // int rowIdn = (i % 64) + 1 ; // for MulRan dataset, Ouster OS1-64 .bin file,  giseop 

            if (rowIdn < 0 || rowIdn >= N_SCAN)
                continue;

            if (rowIdn % downsampleRate != 0)
                continue;

			/* 根据点的水平视角算出点的列号，对应到深度图的横坐标columnIdn */
            float horizonAngle = atan2(thisPoint.x, thisPoint.y) * 180 / M_PI;

            static float ang_res_x = 360.0/float(Horizon_SCAN);
            int columnIdn = -round((horizonAngle-90.0)/ang_res_x) + Horizon_SCAN/2;
            if (columnIdn >= Horizon_SCAN)
                columnIdn -= Horizon_SCAN;

            if (columnIdn < 0 || columnIdn >= Horizon_SCAN)
                continue;

			/* 如果深度图的该位置下已经有数据则略过 */
            if (rangeMat.at<float>(rowIdn, columnIdn) != FLT_MAX)
                continue;

			/* 基于IMU的角速度的积分结果对当前点进行帧内校正 */
            thisPoint = deskewPoint(&thisPoint, laserCloudIn->points[i].time);

            /* 在深度图中填入当前点的深度 */
            rangeMat.at<float>(rowIdn, columnIdn) = range;

			/* 按照行号列号对点云进行重排列，存储到fullCloud中 */
            int index = columnIdn + rowIdn * Horizon_SCAN;
            fullCloud->points[index] = thisPoint;
        }
    }

	/**
	 * 从整理好的点云中提取有效点，顺序一维化，添加到cloudInfo，准备发布。
	 * 在projectPointCloud中点云被按照行列的方式整理成二维方式，其中有些点
	 * 是无效的点，提取过程中这些无效点被滤除。
	 */
    void cloudExtraction()
    {
        int count = 0;
        // extract segmented cloud for lidar odometry
        for (int i = 0; i < N_SCAN; ++i)
        {
            cloudInfo.startRingIndex[i] = count - 1 + 5;

            for (int j = 0; j < Horizon_SCAN; ++j)
            {
                if (rangeMat.at<float>(i,j) != FLT_MAX)
                {
                    // mark the points' column index for marking occlusion later
                    cloudInfo.pointColInd[count] = j;
                    // save range info
                    cloudInfo.pointRange[count] = rangeMat.at<float>(i,j);
                    // save extracted cloud
                    extractedCloud->push_back(fullCloud->points[j + i*Horizon_SCAN]);
                    // size of extracted cloud
                    ++count;
                }
            }
            cloudInfo.endRingIndex[i] = count -1 - 5;
        }
    }
    
	/**
	 * 发布点云cloudInfo
	 */
    void publishClouds()
    {
        cloudInfo.header = cloudHeader;
        cloudInfo.cloud_deskewed  = publishCloud(&pubExtractedCloud, extractedCloud, cloudHeader.stamp, lidarFrame);
        pubLaserCloudInfo.publish(cloudInfo);
    }
};

int main(int argc, char** argv)
{
    ros::init(argc, argv, "lio_sam");

    ImageProjection IP;
    
    ROS_INFO("\033[1;32m----> Image Projection Started.\033[0m");

    ros::MultiThreadedSpinner spinner(3);
    spinner.spin();
    
    return 0;
}
