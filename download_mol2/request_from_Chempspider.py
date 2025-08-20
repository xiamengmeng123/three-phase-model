#爬虫批量下载mol文件
import os
import urllib.request
import requests
import socket
import time

    
for i in range(0,200000):  ##10万开始下
    ##下载2D分子图片
    # url = f'http://www.chemspider.com/ImagesHandler.ashx?id={i}&amp;w=250&amp;h=250'
    #urllib.request.urlretrieve(url,f'{i}.png')
    
    ##批量下载.mol文件
    url = f'http://www.chemspider.com/FilesHandler.ashx?type=str&striph=yes&id={i}'
    file_save_path = 'D:/vspractice/mol_file/200000-250000'    ##保存路径更改

    if not os.path.exists(file_save_path):
        os.mkdir(file_save_path)

    print(f"第{i}个mol文件正在下载......")
    #interval = 2  #设置间隔时间 s
    socket.setdefaulttimeout(30)  #防止超时  单位s
    try:
        urllib.request.urlretrieve(url,f'{file_save_path}/{i}.mol')
        print(f"第{i}个mol文件下载完成！")
        #time.sleep(interval)
    
    except socket.timeout:
        count = 1
        while count <= 5:
            try:
                urllib.request.urlretrieve(url,f'{file_save_path}/{i}.mol')
                print(f"第{i}个mol文件下载完成!")                                                
                break
            except socket.timeout:
                err_info = 'Reloading for %d time'%count if count == 1 else 'Reloading for %d times'%count
                print(err_info)
                count += 1
        if count > 5:
            print("download job failed!")       



    

