# Upload file
curl https://upload.box.com/api/2.0/files/content \
  -H "Authorization: Bearer ACCESS_TOKEN" -X POST \
  -F attributes='{"name":"20200716_065916PM_CDTmolecular_weight_Beniah_10ktest.txt", "parent":{"id":"0"}}' \
  -F file=@20200716_065916PM_CDTmolecular_weight_Beniah_10ktest.txt
# Darn Northwestern's Box settings won't let me generate the access token
