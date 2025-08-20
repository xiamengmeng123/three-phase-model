#!/bin/bash
# 轻量级作业监控与提交守护进程（多分区支持）
# 用法: nohup ./job_daemon.sh <分区名> > /dev/null 2>&1 &

# 配置参数
USER="mmx"
PARTITION="$1"
FEP_DIR="/fs1/home/chengkun_wu/mmx/gromacs"

# 分区特定配置
case $PARTITION in
    cp1) 
        JOB_LIST="$FEP_DIR/job_list_cp1.txt"
        FF="cp1_ff"  # 修改为实际值
        ;;
    cp2)
        JOB_LIST="$FEP_DIR/job_list_cp2.txt"
        FF="cp2_ff"  # 修改为实际值
        ;;
    *)
        echo "无效分区: $PARTITION" >&2
        exit 1
        ;;
esac

# 监控参数
STALL_THRESHOLD=300
SUBMIT_INTERVAL=60
MONITOR_INTERVAL=300

# 轻量级作业列表加载
declare -a JOB_COMMANDS
if [[ -f "$JOB_LIST" ]]; then
    mapfile -t JOB_COMMANDS < "$JOB_LIST"
    TOTAL_JOBS=${#JOB_COMMANDS[@]}
else
    echo "作业列表不存在: $JOB_LIST" >&2
    exit 1
fi

# 初始化提交状态
SUBMITTED=0
LAST_MONITOR_TIME=0

# 主守护进程循环
while true; do
    CURRENT_TIME=$(date +%s)
    
    # === 作业提交逻辑 ===
    if (( SUBMITTED < TOTAL_JOBS )); then
        # 轻量级作业计数
        CURRENT_JOBS=$(yhq -u $USER -p $PARTITION 2>/dev/null | awk 'NR>1' | wc -l)  #过滤掉表头
        AVAILABLE_SLOTS=$((15 - CURRENT_JOBS))
        
        if (( AVAILABLE_SLOTS > 0 )); then
            TO_SUBMIT=$(( AVAILABLE_SLOTS < (TOTAL_JOBS - SUBMITTED) ? AVAILABLE_SLOTS : (TOTAL_JOBS - SUBMITTED) ))
            
            # 批量提交作业
            for ((i = SUBMITTED; i < SUBMITTED + TO_SUBMIT; i++)); do
                eval "${JOB_COMMANDS[$i]}" >/dev/null 2>&1 &
            done
            wait
            
            SUBMITTED=$((SUBMITTED + TO_SUBMIT))
        fi
    fi

    # === 作业监控逻辑 ===
    # 按监控间隔执行
    if (( CURRENT_TIME - LAST_MONITOR_TIME >= MONITOR_INTERVAL )); then
        LAST_MONITOR_TIME=$CURRENT_TIME
        
        # 获取运行中的作业列表
        RUNNING_JOBS=$(yhq -u $USER -h -t R -o "%i" 2>/dev/null)
        
        for JOBID in $RUNNING_JOBS; do
            # 高效获取输出文件路径
            OUTPUT_FILE=$(scontrol show job $JOBID 2>/dev/null | 
                         awk -F= '/StdOut/{print $2; exit}')
            
            if [[ -n "$OUTPUT_FILE" && -f "$OUTPUT_FILE" ]]; then
                LAST_MODIFIED=$(stat -c %Y "$OUTPUT_FILE" 2>/dev/null)
                
                if [[ -n "$LAST_MODIFIED" ]]; then
                    TIME_DIFF=$((CURRENT_TIME - LAST_MODIFIED))
                    
                    # 判断是否卡住
                    if (( TIME_DIFF > STALL_THRESHOLD )); then
                        yhcancel $JOBID >/dev/null 2>&1
                    fi
                fi
            fi
        done
    fi

    # 等待下次检查
    sleep $SUBMIT_INTERVAL
done