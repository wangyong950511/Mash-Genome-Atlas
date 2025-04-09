# 更新脚本到服务器
sudo cp ui.R server.R /srv/shiny-server/myapp/

# 查看日志
ls -lt /var/log/shiny-server/myapp-*.log | head -n 1
tail -f /var/log/shiny-server/myapp-shiny-20250409-135153-39289.log

# 重启shiny-server
sudo systemctl restart shiny-server
