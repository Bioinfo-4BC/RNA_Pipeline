#!/bin/bash
cd {{location}}
samplesCT=({{samplenamesCT}})
samplesST8=({{samplenamesST8}})

bsidsCT=()
bsidsST8=()
other_samples=()
not_fetched=()

for i in "${samplesCT[@]}"; do
  echo "Processing sample: $i"
  attempt_count=0

  while true; do
    bsidCT=$(bs get biosample -n "$i" –terse | grep "Id" | head -1 | grep -Eo '[0-9]{1,}')

    if [ -n "$bsidCT" ]; then
      bsidsCT+=("$bsidCT")
      other_samples+=("$bsidCT")
      break
    else
      echo "No ID found for sample: $i. Retrying after 3 minutes."
      sleep 280
      ((attempt_count++))  # Increment attempt counter

      if [ $attempt_count -eq 5 ]; then
        echo "Moving to next sample for search."
        not_fetched+=("$i")
        break
      fi
    fi

  done
done

for i in "${samplesST8[@]}"; do
  echo "Processing sample: $i"
  attempt_count=0

  while true; do
    bsidST8=$(bs get biosample -n "$i" –terse | grep "Id" | head -1 | grep -Eo '[0-9]{1,}')

    if [ -n "$bsidST8" ]; then
      bsidsST8+=("$bsidST8")
      other_samples+=("$bsidST8")
      break
    else
      echo "No ID found for sample: $i. Retrying after 3 minutes."
      sleep 280
      ((attempt_count++))  # Increment attempt counter

      if [ $attempt_count -eq 5 ]; then
        echo "Moving to next sample for search."
        not_fetched+=("$i")
        break
      fi
    fi

  done
done

printf -v joined '%s,' "${bsidsST8[@]}"
bsidST8=${joined%,}

printf -v joined '%s,' "${bsidsCT[@]}"
bsidCT=${joined%,}

# check if bsidCT is not empty
if [ ${#bsidsCT[@]} -eq 0 ]; then
  echo "No biosample IDs found for CT samples."
else
  echo "Biosample IDs for CT samples: ${bsidsCT[@]}"
  
  {{bscmdct}}
  
fi

# check if bsidST8 is not empty
if [ ${#bsidsST8[@]} -eq 0 ]; then
  echo "No biosample IDs found for ST8 samples."
else
  echo "Biosample IDs for ST8 samples: ${bsidsST8[@]}"
  
  {{bscmdst8}}
  
fi

