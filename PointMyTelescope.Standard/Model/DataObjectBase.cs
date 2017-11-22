using GalaSoft.MvvmLight;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using System.Runtime.InteropServices;

namespace Lunatic.Core
{
   [ComVisible(false)]
   public abstract class DataObjectBase : ObservableObject, INotifyDataErrorInfo
   {
      #region Notify data error
      public event EventHandler<DataErrorsChangedEventArgs> ErrorsChanged;
      private Dictionary<string, List<string>> _errors = new Dictionary<string, List<string>>();

      // get errors by property
      public IEnumerable GetErrors(string propertyName)
      {
         if (!string.IsNullOrWhiteSpace(propertyName)
               && _errors.ContainsKey(propertyName)) {
            return _errors[propertyName];
         }
         return null;
      }

      // has errors
      public bool HasErrors
      {
         get { return (_errors.Count > 0); }
      }


      /// <summary>
      /// Adds the specified error to the errors collection if it is not 
      /// already present, inserting it in the first position if isWarning is 
      /// false. Raises the ErrorsChanged event if the collection changes. 
      /// </summary>
      /// <param name="propertyName"></param>
      /// <param name="error"></param>
      /// <param name="isWarning"></param>
      public void AddError(string propertyName, string error, bool isWarning = false)
      {
         if (!_errors.ContainsKey(propertyName)) {
            _errors[propertyName] = new List<string>();
         }

         if (!_errors[propertyName].Contains(error)) {
            if (isWarning) {
               _errors[propertyName].Add(error);
            }
            else {
               _errors[propertyName].Insert(0, error);
            }
            RaiseErrorsChanged(propertyName);
         }
      }

      /// <summary>
      /// Removes the specified error from the errors collection if it is
      /// present. Raises the ErrorsChanged event if the collection changes.
      /// </summary>
      /// <param name="propertyName"></param>
      /// <param name="error"></param>
      public void RemoveError(string propertyName, string error)
      {
         // remove error
         if (_errors.ContainsKey(propertyName) &&
                _errors[propertyName].Contains(error)) {
            _errors[propertyName].Remove(error);
            if (_errors[propertyName].Count == 0) {
               _errors.Remove(propertyName);
            }
            RaiseErrorsChanged(propertyName);
         }
      }

      /// <summary>
      /// Removes all errors for the given property.
      /// </summary>
      /// <param name="propertyName"></param>
      public void RemoveErrors(string propertyName)
      {
         // remove error
         if (_errors.ContainsKey(propertyName)) {
            _errors.Remove(propertyName);
         }
         RaiseErrorsChanged(propertyName);
      }


      /// <summary>
      /// Raises the ErrorsChanged event.
      /// </summary>
      /// <param name="propertyName"></param>
      public void RaiseErrorsChanged(string propertyName)
      {
         EventHandler<DataErrorsChangedEventArgs> handler = ErrorsChanged;
         // Notify
         if (handler != null) {
            handler(this, new DataErrorsChangedEventArgs(propertyName));
         }
      }
      #endregion

   }
}
