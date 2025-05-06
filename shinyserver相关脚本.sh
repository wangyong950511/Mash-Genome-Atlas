# 更新脚本到服务器
sudo cp ui.R server.R global.R /srv/shiny-server/myapp/

# 查看日志
tail -f $(ls -t /var/log/shiny-server/myapp-*.log | head -n 1)

# 重启shiny-server
sudo systemctl restart shiny-server

## 网站维护
# 防火墙状态
sudo tail -n 100 /var/log/nginx/access.log
