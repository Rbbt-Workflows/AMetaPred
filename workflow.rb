require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

Workflow.require_workflow "Sequence"
Workflow.require_workflow "DbNSFP"

module PMut
  extend Workflow

  input :mutations, :array, "Mutated isoforms"
  task :predict => :tsv do |mutations|
    index = Organism.protein_identifiers("Hsa/feb2014").index :target => "UniProt/SwissProt Accession", :persist => true
    uniprot_mutations = mutations.collect do |mi|
      protein, change = mi.split(":")
      uniprot = index[protein]
      next if uniprot.nil?
      [uniprot, change, mi]
    end.compact.uniq

    acc = nil
    TSV.traverse uniprot_mutations, :type => :array, :bar => self.progress_bar("PMut REST queries") do |uniprot,change,mi|
      ref, pos, alt = change.match(/(.)(\d+)(.+)/).values_at 1, 2, 3
      url = "https://mmb.irbbarcelona.org/PMut/uniprot/#{uniprot}/#{pos}/#{alt}/features.csv"
      begin
        tsv = TSV.open(url, :sep => ',', :header_hash => '', :type => :list)
      rescue
        next
      end
      tsv.add_field "Mutation" do 
        mi
      end
      tsv = tsv.reorder "Mutation"
      if acc.nil?
        acc = tsv
      else
        acc.merge!(tsv)
      end
    end
    
    acc
  end
end

module AMetaPred
  extend Workflow

  task :parse_mutations => :tsv do
    input = Rbbt.data["DDD_input.csv"].tsv :type => :list, :header_hash => '', :sep => ','
    input.fields = input.fields.collect{|f| f.sub('#', '') }
    input.key_field = "Num"
    output = Rbbt.data["DDD_output.csv"].tsv :type => :list, :header_hash => '', :sep => ','
    output.key_field = "Num"
    input.add_field "Genomic Mutation" do |n,values|
      m = values["Variant"]
      chr, pos, ref, alt = m.match(/(.+)\(.+:g\.(\d+)(.+)>(.+)/).values_at 1, 2, 3, 4
      chr.sub!('Chr', '')
      [chr, pos, alt] * ":"
    end
    
    input.add_field "Fujitsu" do
    end
    input.attach output
  end

  dep :parse_mutations
  task :genomic_mutations => :array do
    step(:parse_mutations).load.column("Genomic Mutation").values.flatten
  end

  dep :genomic_mutations
  dep_task :mutated_isoforms, Sequence, :mutated_isoforms_fast, :mutations => :genomic_mutations, :organism => "Hsa/feb2014", :principal => false

  dep :mutated_isoforms
  task :mis => :array do 
    step(:mutated_isoforms).load.values.flatten.uniq
  end

  dep :mis
  dep_task :mi_predictions, DbNSFP, :predict, :mutations => :mis

  dep :mis
  dep_task :pmut, PMut, :predict, :mutations => :mis

  dep :pmut
  dep :mi_predictions
  task :add_pmut => :tsv do
    tsv = step(:mi_predictions).load
    pmut = step(:pmut).load
    pmut.add_field "PMut" do |k,v|
      v["pred_disease"].downcase == "true" ? "D" : "B"
    end
    pmut = pmut.slice(["PMut"])
    pmut.key_field = tsv.key_field
    tsv.attach pmut
  end


  dep :mi_predictions
  task :values => :tsv do
    tsv = step(:mi_predictions).load
    tsv.fields.each do |f|
      iii [f, Misc.counts(tsv.column(f).values.flatten)]
    end
  end

  dep :add_pmut
  dep :mutated_isoforms
  dep :parse_mutations
  task :result => :tsv do
    tsv = step(:parse_mutations).load
    predictions = step(:add_pmut).load
    mis = step(:mutated_isoforms).load

    predictions.fields.each do |predictor|
      predictor_values = predictions.column(predictor)
      tsv.add_field predictor do |key,values|
        mutation = values["Genomic Mutation"]
        mutation_mis = mis[mutation]
        raise mutation if mutation_mis.nil?
        preds = predictor_values.values_at *mutation_mis

        pathogenic = false
        pathogenic = true if preds.include?("D")
        pathogenic = true if preds.include?("H")
        pathogenic = true if preds.include?("M")
        pathogenic
      end
    end
    tsv
  end

  dep :result
  input :threshold, :float, "Score threshold for Fujitsu", 0.6
  input :revel_threshold, :float, "Score threshold for REVEL", 0.5
  input :clinpred_threshold, :float, "Score threshold for ClinPred", 0.5
  task :add_fujitsu => :tsv do |threshold,revel_threshold,clinpred_threshold|
    tsv = step(:result).load
    tsv.add_field "REVEL" do |k,v|
      v["REVEL_SCORE"].to_f > revel_threshold
    end
    tsv.add_field "ClinPred_ORIG" do |k,v|
      v["CLINPRED_SCORE"].to_f > clinpred_threshold
    end
    tsv.add_field "Fujitsu-XAIscoring" do |k,v|
      v["XAIscoring"].to_f > threshold
    end
  end

  dep :add_fujitsu
  task :evaluate => :tsv do
    result = step(:add_fujitsu).load
    predictors = result.fields[7..-1]
    tsv = TSV.setup({}, "Predictor~TP,TN,FP,FN,Precision,Recall,F-score,MCC#:type=:list#:cast=:to_f")

    predictors.each do |predictor|
      next if predictor == "Aloft_pred"
      tp = 0
      tn = 0
      fp = 0
      fn = 0
      result.through do |mut,values|
        label = values["PATHOGENICITY"] == 'PATHOGENIC'
        pred_value = values[predictor]
        next if pred_value.nil? || pred_value.empty?
        pred = pred_value.to_s == 'true'
        if label
          if pred
            tp += 1
          else
            fn += 1
          end
        else
          if pred
            fp += 1
          else
            tn += 1
          end
        end
      end
      recall = tp.to_f / (tp + fn)
      precision = tp.to_f / (tp + fp)

      f_score = 2 * precision * recall / (precision + recall)
      
      mcc = ((tp * tn) - (fp*fn)) / Math.sqrt( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

      tsv[predictor] = [tp, tn, fp, fn, precision, recall,f_score,mcc]
    end
    tsv
  end

  dep :evaluate
  task :json => :json do
    orig = step(:evaluate).load
    tsv = orig.annotate({})
    
    orig.each do |k,v|
     tsv[k.sub("_pred",'_P')] = v
    end

    community_id = "AMetaPred"
    event_id = "AMetaPred-DDD"
    challenge_id = "RD"

    date = Time.now
    participants = tsv.collect do |participant_id,values|
      {
        "_id" => [community_id, participant_id] * ":",
        "challenge_id" => [
          challenge_id
        ],
        "community_id" => community_id,
        "datalink" => {
          "attrs" => [
            "archive"
          ],
          "status" => "ok",
          "validation_date" => date
        },
        "participant_id" => participant_id,
        "type" => "participant"
      }
    end

    assessments = [] 
    tsv.fields.each do |metric_id|
      assessments += tsv.column(metric_id).collect do |participant_id,value|
        {
          "_id" => "#{community_id}:#{[challenge_id, participant_id, metric_id, "A"] * "_"}",
          "challenge_id" => challenge_id,
          "community_id" => community_id,
          "metrics" => {
            "metric_id" => metric_id,
            "value" => value
          },
          "participant_id" => participant_id,
          "type" => "assessment"
        }
      end
    end

    x_field = "Precision"
    y_field = "Recall"
    aggregation_data = tsv.keys.collect do |participant_id|
      {
        "metric_x" => 218462,
        "metric_y" => 191952,
        "participant_id" => participant_id
      }
    end

    aggregation = {
      "_id" => "#{community_id}:#{[challenge_id, 'agg', [x_field, y_field] * "+"] * "_"}",
        "challenge_ids" => [
          challenge_id
        ],
        "community_id" => community_id,
          "datalink" => {
            "inline_data" => {
              "challenge_participants" => aggregation_data,
              "visualization" => {
                "type" => "2D-plot",
                "x_axis" => x_field,
                "y_axis" => y_field
              }
            }
          },
          "type" => "aggregation"
    }

    challenge = {
      "id" => challenge_id,
      "participants" => tsv.keys
    }
    participants + assessments + [aggregation] + [challenge]
  end

  dep :evaluate
  task :tool_jsons => :array do
    orig = step(:evaluate).load
    tsv = orig.annotate({})
    
    community_id = "AMetaPred"
    event_id = "AMetaPred-DDD"
    challenge_id = "RD"

    orig.each do |k,v|
     tsv[k.sub("_pred",'')] = v
    end

    json_dir = file('jsons')
    tsv.keys.each do |participant_id|
      tool_access_type = 'command-line'
      link = "https://sites.google.com/site/jpopgen/dbNSFP"
      contact = "Miguel.Vazquez.Garcia"
      info = {
        "_id" => participant_id,
        "_schema" => "https://www.elixir-europe.org/excelerate/WP2/json-schemas/1.0/Tool",
        "community_ids" => [community_id],
        "name" => participant_id,
        "description" => "Tool #{participant_id} from DbNSFP",
        "is_automated" => false,
        "tool_contact_ids" => [
          contact
        ],
        "status" => "online",
        "references" => [  ],
        "tool_access" => [
          {
            "tool_access_type" => tool_access_type,
            "link" => link,
            "techdocs" => [{"uri":"https://github.com/MathCancer/PhysiCell/tree/master/documentation","label":"source"}]
          }
        ]

      }  
      Open.write(json_dir["Tool.#{participant_id}.json"], info.to_json)
    end

    json_dir.glob("*.json")
  end

  dep :evaluate
  task :metric_jsons => :array do
    orig = step(:evaluate).load
    tsv = orig.annotate({})
    
    community_id = "AMetaPred"
    event_id = "AMetaPred-DDD"
    challenge_id = "RD"

    orig.each do |k,v|
     tsv[k.sub("_pred",'')] = v
    end

    json_dir = file('jsons')
    tsv.fields.each do |metric_id|
      info = {"_id" => metric_id,
       "_schema" => "https://www.elixir-europe.org/excelerate/WP2/json-schemas/1.0/Metrics",
       "title": metric_id,
       "metrics_contact_ids": ["Miguel.Vazquez.Garcia"],
      }
      Open.write(json_dir["Tool.#{metric_id}.json"], info.to_json)
    end

    json_dir.glob("*.json")
  end


  dep :evaluate
  task :dataset_json => :array do
    orig = step(:evaluate).load
    tsv = orig.annotate({})
    
    community_id = "AMetaPred"
    event_id = "AMetaPred-DDD"
    challenge_id = "RD"

    orig.each do |k,v|
     tsv[k.sub("_pred",'')] = v
    end

    json_dir = file('jsons')
    tsv.fields.each do |metric_id|
      {"_id" => metric_id,
       "_schema" => "https://www.elixir-europe.org/excelerate/WP2/json-schemas/1.0/Metrics",
       "title": metric_id,
       "metrics_contact_ids": ["Miguel.Vazquez.Garcia"],
      }

    end

    json_dir.glob("*.json")
  end

end

#require 'MODULE/tasks/basic.rb'

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

