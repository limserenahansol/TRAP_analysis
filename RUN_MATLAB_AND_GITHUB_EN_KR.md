# How to run in MATLAB / MATLAB에서 실행하기

## English

1. Open **MATLAB**.
2. Set the current folder to the pipeline root (folder that contains `trap_config.m` and `init_TRAP_pipeline.m`), e.g.  
   `C:\Users\hsollim\Desktop\TRAP_pipeline`
3. In the Command Window:

```matlab
init_TRAP_pipeline
RUN_PIPELINE_ALL
```

4. **First-time setup:** edit `trap_config.m` (paths), create `TRAP_cohort_CSVs.txt` and `TRAP_sample_manifest.csv` as in `README.md`.
5. **Run only Steps 6–9 (after 1–5 already done):**

```matlab
init_TRAP_pipeline
trap_run_phase_AP_contrasts      % Step 6
trap_run_phase_delta_screening    % Step 6b
trap_run_directional_AP_scenario_folders  % Step 7
trap_run_phase_delta_within_group         % Step 8
trap_run_step9_forebrain_exclude_bs_cb    % Step 9
```

Outputs go under **`TRAP_OUTPUT/`**.

---

## 한국어

1. **MATLAB** 실행.
2. 현재 폴더를 `trap_config.m`이 있는 **TRAP_pipeline 루트**로 설정 (예: `C:\Users\hsollim\Desktop\TRAP_pipeline`).
3. 명령 창에서:

```matlab
init_TRAP_pipeline
RUN_PIPELINE_ALL
```

4. **처음 한 번:** `trap_config.m`, `TRAP_cohort_CSVs.txt`, `TRAP_sample_manifest.csv` 설정 (`README.md` 참고).
5. **6~9단계만** (1~5 이미 돌린 뒤) 위 코드 블록처럼 각 `trap_run_*` 호출.

결과는 **`TRAP_OUTPUT/`** 아래에 저장됩니다.

---

# Push to GitHub / GitHub에 올리기

## English

1. Install **Git** ([git-scm.com](https://git-scm.com/)). On GitHub, create a **new repository** (empty, no README required).

2. In **PowerShell** or **Git Bash**, go to your pipeline folder:

```powershell
cd C:\Users\hsollim\Desktop\TRAP_pipeline
```

3. If this folder is **not** yet a git repo:

```powershell
git init
git add .
git commit -m "Initial TRAP pipeline"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO.git
git push -u origin main
```

Replace `YOUR_USERNAME/YOUR_REPO` with your GitHub user and repo name.

4. **Large CSV / TRAP_OUTPUT:** add to `.gitignore` if you do not want them online:

```
TRAP_OUTPUT/
*.csv
```

Commit `.gitignore`, then `git add` / `git commit` / `git push` again.

5. **Updates later:**

```powershell
git add .
git commit -m "Describe your change"
git push
```

---

## 한국어

1. **Git** 설치 후 GitHub에서 **새 저장소** 생성.

2. PowerShell에서 파이프라인 폴더로 이동:

```powershell
cd C:\Users\hsollim\Desktop\TRAP_pipeline
```

3. 아직 git이 없으면:

```powershell
git init
git add .
git commit -m "Initial TRAP pipeline"
git branch -M main
git remote add origin https://github.com/본인아이디/저장소이름.git
git push -u origin main
```

4. **용량 큰 CSV·결과 폴더**는 올리지 않으려면 `.gitignore`에 `TRAP_OUTPUT/`, `*.csv` 등 추가 후 커밋.

5. **이후 수정** 후에는 `git add .` → `git commit -m "메시지"` → `git push`.

---

**Personal access token:** GitHub는 비밀번호 대신 **PAT**로 `git push`할 수 있습니다. GitHub → Settings → Developer settings → Personal access tokens.
