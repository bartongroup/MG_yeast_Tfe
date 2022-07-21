rsync -rvm --include='*/' --include='doc/analysis.html' --include='doc/image/*' --include='tab/*' \
  --exclude='*' . cluster:/cluster/gjb_lab/mgierlinski/public_html/yeast_tfe
